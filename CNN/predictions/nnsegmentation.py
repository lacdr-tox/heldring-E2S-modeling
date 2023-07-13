print("Initializating...")

import sys
import os
from os import listdir, mkdir, path, walk
from os.path import isfile, isdir, join, basename
import math
import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
from torch.optim.lr_scheduler import ReduceLROnPlateau
import torchvision
import torchvision.transforms as transforms
import numpy as np
import matplotlib.pyplot as plt
from zipfile import ZipFile
import cv2
from PIL import Image
from torch.autograd import Variable
from scipy.ndimage.interpolation import map_coordinates
from scipy.ndimage.filters import gaussian_filter
from scipy.ndimage import distance_transform_edt, grey_closing
from scipy.ndimage.filters import gaussian_filter
from skimage.segmentation import watershed
from skimage.measure import label
from skimage import util
from sklearn.utils import shuffle
from tqdm import tqdm
import tifffile as tiff
from tifffile import imsave
import time

print("All libraries are successfully loaded!")

print("This is the script: %s" %(sys.argv[0]))
print("This is the path to the models: %s" %(sys.argv[1]))
print("This is the input directory: %s" %(sys.argv[2]))
print("This is the output directory: %s" %(sys.argv[3]))


#Normalize
def normalize(img, low, high):
  min=img.min()
  max=img.max()
  a=img-min
  b=max-min
  return (high-low)*(a/b)+low

def crop(image):
    width, height = image.size 
    left = 0
    top = 0
    right = width-(width%16)
    bottom = height-(height%16)
    image = image.crop((left, top, right, bottom))
    return image

def process(image, norm):
    image.convert('I;16')
    image = crop(image)
    image = np.asarray(image)
    if norm:
      image = normalize(image,-1,1)
    w=image.shape[0]
    h=image.shape[1]
    image = np.reshape(image, (1, w, h))
    return image,w,h

def elastic_transform(image, alpha, sigma, random_state=None):
    """Elastic deformation of images as described in [Simard2003]_.
    .. [Simard2003] Simard, Steinkraus and Platt, "Best Practices for
       Convolutional Neural Networks applied to Visual Document Analysis", in
       Proc. of the International Conference on Document Analysis and
       Recognition, 2003.
    """
    assert len(image.size)==2

    if random_state is None:
        random_state = np.random.RandomState(None)

    shape = image.size

    dx = gaussian_filter((random_state.rand(*shape) * 2 - 1), sigma, mode="constant", cval=0) * alpha
    dy = gaussian_filter((random_state.rand(*shape) * 2 - 1), sigma, mode="constant", cval=0) * alpha

    x, y = np.meshgrid(np.arange(shape[0]), np.arange(shape[1]), indexing='ij')
    indices = np.reshape(x+dx, (-1, 1)), np.reshape(y+dy, (-1, 1))
    
    return map_coordinates(image, indices, order=1).reshape(shape)

#rotate, mirror, elastic deformation
def augment(image, original):
  out=[]
  out.append(image.rotate(90))
  out.append(image.rotate(180))
  out.append(image.rotate(270))
  mir=image.transpose(Image.FLIP_LEFT_RIGHT)
  out.append(mir)
  out.append(mir.rotate(90))
  out.append(mir.rotate(180))
  out.append(mir.rotate(270))

  rndnr=7
  transformed_image=Image.fromarray(elastic_transform(image,40,6,random_state=np.random.RandomState(rndnr)))
  out.append(transformed_image.rotate(90))
  out.append(transformed_image.rotate(180))
  out.append(transformed_image.rotate(270))
  transformed_image_mir=transformed_image.transpose(Image.FLIP_LEFT_RIGHT)
  out.append(transformed_image_mir)
  out.append(transformed_image_mir.rotate(90))
  out.append(transformed_image_mir.rotate(180))
  out.append(transformed_image_mir.rotate(270))

  rndnr=8
  new_signal=elastic_transform(image,40,6,random_state=np.random.RandomState(rndnr))
  new_signal=Image.fromarray(new_signal)
  out.append(new_signal.rotate(90))
  out.append(new_signal.rotate(180))
  out.append(new_signal.rotate(270))
  new_signal_mir=new_signal.transpose(Image.FLIP_LEFT_RIGHT)
  out.append(new_signal_mir)
  out.append(new_signal_mir.rotate(90))
  out.append(new_signal_mir.rotate(180))
  out.append(new_signal_mir.rotate(270))
  return out


def load_training(input_dir,cdp_dir,ndp_dir):
    print("Loading images from: "+input_dir)
    start_time = time.time()
    
    inputf = [Image.open(join(input_dir, f)) for f in sorted(listdir(input_dir)) if isfile(join(input_dir, f))]
    inputc = [Image.open(join(cdp_dir, f)) for f in sorted(listdir(cdp_dir)) if isfile(join(cdp_dir, f))]
    inputn = [Image.open(join(ndp_dir, f)) for f in sorted(listdir(ndp_dir)) if isfile(join(ndp_dir, f))]
    
    print("amount of files: "+str(len(inputf)))
    input=[]
    length=len(inputf)
    #augment
    for i in range(length):      
      inputf.extend(augment(inputf[i],True))
      inputc.extend(augment(inputc[i],False))
      inputn.extend(augment(inputn[i],False))
    print("amount of training after augmentation: "+str(len(inputf)))
    for i in range(len(inputf)):
      inputf[i],w,h=process(inputf[i],True)
      inputc[i],w,h=process(inputc[i],True)
      inputn[i],w,h=process(inputn[i],True)

    inputf, inputc, inputn = shuffle(np.array(inputf), np.array(inputc), np.array(inputn))
    train_in=torch.from_numpy(inputf).float()
    train_cdp=torch.from_numpy(inputc).float()
    train_ndp=torch.from_numpy(inputn).float()
    print("--- %s seconds ---" % (time.time() - start_time))
    print("Images loaded")
    return train_in, train_cdp, train_ndp

def save_img(array,name):
  array=normalize(array,0,65535)
  tiff.imsave(name,array)

class UNet(nn.Module):
    def __init__(
        self,
        decoder,
        in_channels=1,
        n_classes=1,
        depth=5,
        wf=6,
        padding=True,
        batch_norm=True,
        up_mode='upconv',
    ):
        super(UNet, self).__init__()        
        assert up_mode in ('upconv', 'upsample')
        self.decoder=decoder
        assert decoder in ('ndp', 'cdp')
        self.padding = padding
        self.depth = depth
        prev_channels = in_channels
        
        self.encoderConv = nn.ModuleList()
        self.pooling = nn.ModuleList()
        for i in range(depth):
            self.encoderConv.append(
                UNetConvBlock(prev_channels, 2 ** (wf + i), padding, batch_norm)
            )
            if i+1<=depth-1:
                self.pooling.append(
                    UNetDownBlock(2 ** (wf + i), 2 ** (wf + i), padding, batch_norm)
                )          
            prev_channels = 2 ** (wf + i)
        
        if (decoder=="ndp"):
          self.decoder1Conv = nn.ModuleList()
          self.decoder1Upconv = nn.ModuleList()          
          for i in reversed(range(depth - 1)):
              self.decoder1Upconv.append(
                  UNetUpBlock(prev_channels, 2 ** (wf + i), up_mode, padding, batch_norm)
              )
              self.decoder1Conv.append(
                  UNetConvBlock(prev_channels, 2 ** (wf + i), padding, batch_norm)
              )
              prev_channels = 2 ** (wf + i)          
          self.decoder1Conv.append(nn.Conv2d(prev_channels, n_classes, kernel_size=1))

        else:          
          self.decoder2Conv = nn.ModuleList()
          self.decoder2Upconv = nn.ModuleList()
          for i in reversed(range(depth - 1)):              
              self.decoder2Upconv.append(
                  UNetUpBlock(prev_channels, 2 ** (wf + i), up_mode, padding, batch_norm)
              )            
              self.decoder2Conv.append(
                  UNetConvBlock(prev_channels, 2 ** (wf + i), padding, batch_norm)
              )            
              prev_channels = 2 ** (wf + i)
          self.decoder2Conv.append(nn.Conv2d(prev_channels, n_classes, kernel_size=1))
        

    def forward(self, x):
        blocks = []
        for i, down in enumerate(self.encoderConv):
            x = down(x)
            if i != len(self.encoderConv) - 1:
                pool=self.pooling[i]
                blocks.append(x)           
                x=pool(x)

        if (self.decoder=="ndp"):
          for i, up in enumerate(self.decoder1Upconv):
              conv=self.decoder1Conv[i]
              x = up(x, blocks[-i - 1])
              x=conv(x)
          last=self.decoder1Conv[self.depth-1]
        else:
          for i, up in enumerate(self.decoder2Upconv):
              conv=self.decoder2Conv[i]
              x = up(x, blocks[-i - 1])
              x=conv(x)
          last=self.decoder2Conv[self.depth-1]
        return last(x)


class UNetConvBlock(nn.Module):
    def __init__(self, in_size, out_size, padding, batch_norm):
        super(UNetConvBlock, self).__init__()
        conv = []

        conv.append(nn.Conv2d(in_size, out_size, kernel_size=3, padding=int(padding)))
        conv.append(nn.ReLU())
        if batch_norm:
            conv.append(nn.BatchNorm2d(out_size))

        conv.append(nn.Conv2d(out_size, out_size, kernel_size=3, padding=int(padding)))
        conv.append(nn.ReLU())
        if batch_norm:
            conv.append(nn.BatchNorm2d(out_size))

        self.conv = nn.Sequential(*conv)

    def forward(self, x):
        out = self.conv(x)
        return out

class UNetDownBlock(nn.Module):
    def __init__(self, in_size, out_size, padding, batch_norm):
        super(UNetDownBlock, self).__init__()
        conv_pool = []

        conv_pool.append(nn.Conv2d(out_size, out_size, kernel_size=3, stride=2, padding=int(padding)))
        conv_pool.append(nn.ReLU())
        if batch_norm:
            conv_pool.append(nn.BatchNorm2d(out_size))

        
        self.conv_pool = nn.Sequential(*conv_pool)

    def forward(self, x):
        out = self.conv_pool(x)
        return out


class UNetUpBlock(nn.Module):
    def __init__(self, in_size, out_size, up_mode, padding, batch_norm):
        super(UNetUpBlock, self).__init__()
        up = []
        if up_mode == 'upconv':
            up.append(nn.ConvTranspose2d(in_size, out_size, kernel_size=2, stride=2))
            
        elif up_mode == 'upsample':
            up.append(nn.Upsample(mode='bilinear', scale_factor=2))
            up.append(nn.Conv2d(in_size, out_size, kernel_size=1))
        self.up=nn.Sequential(*up)
        self.norm=nn.BatchNorm2d(out_size)

    def center_crop(self, layer, target_size):
        _, _, layer_height, layer_width = layer.size()
        diff_y = (layer_height - target_size[0]) // 2
        diff_x = (layer_width - target_size[1]) // 2
        return layer[
            :, :, diff_y : (diff_y + target_size[0]), diff_x : (diff_x + target_size[1])
        ]

    def forward(self, x, bridge):
        up = self.up(x)
        norm=self.norm(up)
        out=torch.cat([norm,bridge], 1)

        return out

def train(model, output_name, train_in, train_out,epochs,batch_size,learning_rate,val_percent,optimizer,scheduler):
  criterion = nn.SmoothL1Loss()
  best_val_loss=1000
  n_val = int(len(train_in) * val_percent)
  n_train = len(train_in) - n_val

  val_in=train_in[n_train:]
  val_out=train_out[n_train:]
  
  train_in=train_in[:n_train]
  train_out=train_out[:n_train]

  for epoch in range(1, epochs+1):
      
      print('Epoch:', epoch)
      model.train(True)  # Set model to training mode    
      permutation = torch.randperm(train_in.size()[0])
      training_loss = []      
      for i in tqdm(range(0,train_in.size()[0], batch_size)):
          indices = permutation[i:i+batch_size]
          batch_x, batch_y = train_in[indices], train_out[indices]
          batch_x, batch_y = batch_x.cuda(),batch_y.cuda()
          output = model(batch_x)
          loss = criterion(output,batch_y)
          training_loss.append(loss.item())
          optimizer.zero_grad()
          loss.backward()
          optimizer.step()
          
      training_loss = np.average(training_loss)      
      model.train(False)  # Set model to evaluate mode    
      permutation = torch.randperm(val_in.size()[0])
      val_loss = []
      for i in range(0,val_in.size()[0], batch_size):
          indices = permutation[i:i+batch_size]
          batch_x, batch_y= val_in[indices], val_out[indices]  
          batch_x, batch_y = batch_x.cuda(),batch_y.cuda()
          output = model(batch_x)
          loss = criterion(output,batch_y)
          val_loss.append(loss.item())
      val_loss = np.average(val_loss)      
      scheduler.step(val_loss)
      print('\n training loss: \t', training_loss, '\t val loss: \t', val_loss)      
      if    val_loss < best_val_loss:
        best_val_loss=val_loss
        print("saving model")
        torch.save(model.state_dict(), output_name)

  print('Finished Training: ', output_name)

def load_images(path,batch_size):
    print("Loading images from: "+path)
    start_time = time.time()
    files = [join(path, f) for f in sorted(listdir(path)) if isfile(join(path, f))]
    print("amount of files: "+str(len(files)))
    images=[]
    for i in range(math.ceil(len(files)/batch_size)):
      b=[]
      for j in range(batch_size):
        iter=batch_size*i+j
        if iter<len(files):
          f=files[iter]
          image = Image.open(f)
          image,w,h=process(image,True) 
          image=np.reshape(image, (1, w, h))
          b.append(image)
      images.append(b)
    print("--- %s seconds ---" % (time.time() - start_time))
    print("Images loaded")
    return images

def save_pred(pred,name):
  start_time = time.time()
  print("saving predicted imgs")
  preds=[]
  for i,p in enumerate(pred):
    p=p[0]
    p=normalize(p,0,1)
    p = Image.fromarray(p)   
    preds.append(p)
  print("--- %s seconds ---" % (time.time() - start_time))
  print("saving done")
  return preds

def remap( arr, ):
  max_id       = np.max( arr )
  mapping      = np.random.randint( low=1, high=255, size=(int(max_id+1)), dtype=np.uint8)
  mapping[0]=0
  mapped_array = mapping[arr]
  return mapped_array

def example(i):
  labels = prediction[0]
  names = prediction[1]
  images = prediction[2]
  cdp_smooth = prediction[3]
  ndp_square = prediction[4]
  binmask = prediction[5]
  seeds = prediction[6]
  label=remap(labels[i])
  print(names[i])
  fig, axes = plt.subplots(nrows=2,ncols=3, figsize=(18, 12), sharex=True, sharey=True)
  ax = axes.ravel()

  ax[0].imshow(images[i], cmap=plt.cm.gray)
  ax[0].set_title('Original')
  ax[1].imshow(cdp_smooth[i], cmap=plt.cm.gray)
  ax[1].set_title('Smoothed cdp')
  ax[2].imshow(ndp_square[i], cmap=plt.cm.gray)
  ax[2].set_title('Squared ndp')
  ax[3].imshow(binmask[i], cmap=plt.cm.gray)
  ax[3].set_title('Binary mask')
  ax[4].imshow(seeds[i], cmap=plt.cm.gray)
  ax[4].set_title('Seeds')
  ax[5].imshow(label, cmap='nipy_spectral', interpolation='nearest')
  ax[5].set_title('Watershed')

  for a in ax:
      a.set_axis_off()
      a.set_aspect('equal')

  fig.tight_layout()
  plt.subplots_adjust(wspace=0.1, hspace=0.1)
  plt.show()
  print("Dino tries!")
  
def load_model(model,PATH):
  print("loading model: ",PATH)
  device=torch.device('cuda')
  #pretrained_dict = torch.load(PATH, map_location=torch.device('cpu'))
  pretrained_dict = torch.load(PATH, map_location=torch.device('cpu'))
  print("pretrained model loaded")
  model_dict = model.state_dict()
  # 1. filter out unnecessary keys
  pretrained_dict = {k: v for k, v in pretrained_dict.items() if k in model_dict}
  # 2. overwrite entries in the existing state dict
  model_dict.update(pretrained_dict) 
  # 3. load the new state dict
  model.load_state_dict(model_dict)
  model.to(device)
  print("model loaded")
  return model

def to_uint(img):  
  return cv2.normalize(img, img, 0, 65535, cv2.NORM_MINMAX).astype("uint16")

def to_int(img):
  return cv2.normalize(img, img, 0, 255, cv2.NORM_MINMAX).astype("uint8")

def add_imgs(imgs):
  res=imgs[0]
  for i in range(len(imgs)-1):
    cv2.add(res, imgs[i+1], res)
  return res

def preprocess(seg_dir,cdp_dir,ndp_dir):
  print("preprocessing...")
  if not path.exists(cdp_dir):
    mkdir(cdp_dir)
    mkdir(ndp_dir)
  files=[join(seg_dir, f) for f in sorted(listdir(seg_dir)) if isfile(join(seg_dir, f))]
  print("amount of files: ",len(files))
  for iter,f in enumerate(files):
    wshed = cv2.imread(f, False)
    
    cells=[(wshed==i) for i in np.delete(np.unique(wshed), 0, 0)]
    cdps=[]
    for i,cell in enumerate(cells):
      dt=distance_transform_edt(cell)  
      normalizedImg = cv2.normalize(dt, dt, 0, 1, cv2.NORM_MINMAX)      
      cdps.append(normalizedImg)  
    cdp=add_imgs(cdps)
    name="cdp"+('%02d' % iter)+".tif"
    imsave(join(cdp_dir, name),to_uint(cdp))

    ndps=[]
    name="ndp"+('%02d' % iter)+".tif"
    for cell in cells:
      cell=np.array(cell,dtype="uint8")
      retval,cell=cv2.threshold(cell,0,255,cv2.THRESH_BINARY)
      retval,wshed=cv2.threshold(wshed,0,255,cv2.THRESH_BINARY)      
      res=wshed-cell
      res=util.invert(res)    
      res=distance_transform_edt(res)      
      mask=util.invert(cell)
      mask=mask//255  
      mask=(mask==1)
      res[mask] = 0
      res = cv2.normalize(res, res, 0, 1, cv2.NORM_MINMAX)      
      res=util.invert(res)
      res[mask] = -1
      res=res+1
      ndps.append(res)
    ndp=add_imgs(ndps)
    ndp=grey_closing(ndp,(3,3))
    ndp=np.power(ndp, 10)
    imsave(join(ndp_dir, name),to_uint(ndp))
  print("preprocessing done")

def train_models(input_dir,seg_dir, epochs, batch_size):
  cdp_dir="/content/cdp/"
  ndp_dir="/content/ndp/"
  preprocess(seg_dir,cdp_dir,ndp_dir)
  train_in, train_cdp, train_ndp = load_training(input_dir,cdp_dir,ndp_dir)
  assert train_in.size()==train_cdp.size()  
  learning_rate = 8e-4 #8e-4 original
  val_percent=0.1  

  optimizer = optim.Adam(ndp_model.parameters(), lr=learning_rate)
  scheduler = ReduceLROnPlateau(optimizer, mode='min', factor=0.25, patience=10, threshold=0.0001, threshold_mode='rel', cooldown=0, min_lr=8e-5, eps=1e-08, verbose=True)
  train(ndp_model, "ndp_model.pt",train_in, train_ndp, epochs,batch_size,learning_rate,val_percent,optimizer,scheduler)
  optimizer = optim.Adam(cdp_model.parameters(), lr=learning_rate)
  scheduler = ReduceLROnPlateau(optimizer, mode='min', factor=0.25, patience=10, threshold=0.0001, threshold_mode='rel', cooldown=0, min_lr=8e-5, eps=1e-08, verbose=True)
  train(cdp_model, "cdp_model.pt",train_in, train_cdp, epochs,batch_size,learning_rate,val_percent,optimizer,scheduler)

def predict(input_dir, output_dir, tmask, tmarker):
  device=torch.device('cuda')
  path=input_dir
  save_path=output_dir
  if not os.path.exists(save_path):
    os.mkdir(save_path)
  batch_size=10
  images = load_images(path,batch_size)

  print("starting prediction")
  start_time = time.time()
  ndps=[]
  cdps=[]
  ndp_model.train(False)  # Set model to evaluate mode
  cdp_model.train(False)  # Set model to evaluate mode
  with torch.no_grad():
    for i,chunk in enumerate(images):
      torch.cuda.empty_cache()
      print("predicting chunk "+str(i))
      chunk = np.array(chunk)    
      chunk=torch.from_numpy(chunk).float()
      chunk = Variable(chunk, requires_grad=False)
      chunk = chunk.to(device)    
      ndp = ndp_model(chunk)
      cdp=cdp_model(chunk)
      ndps.append(ndp.to('cpu'))
      cdps.append(cdp.to('cpu'))
  print("--- %s seconds ---" % (time.time() - start_time))
  print("prediction done")

  ndp=np.concatenate(ndps)
  cdp=np.concatenate(cdps)
  ndp=save_pred(ndp,"ndp_raw")
  cdp=save_pred(cdp,"cdp_raw")

  sigma=1.5  
  ndp_smooth=[]
  ndp_square=[]
  cdp_smooth=[]
  binmask=[]
  seeds=[]
  labels=[]
  for img in ndp:
    filtered=gaussian_filter(img,sigma)
    ndp_smooth.append(filtered)
    ndp_square.append(np.square(filtered))
  for i,img in enumerate(cdp):
    squared=ndp_square[i]    
    filtered=gaussian_filter(img,sigma)
    cdp_smooth.append(filtered)
    mask=filtered>tmask
    binmask.append(mask)
    seed=label((filtered-squared)>tmarker)
    seeds.append(seed)
    wshed=watershed(-filtered, markers=seed, mask=mask, watershed_line=True)  
    labels.append(wshed) 

  names=[]
  files=[]
  for f in sorted(listdir(path)):
    if isfile(join(path, f)):
      files.append(join(path, f))
      names.append(f)
  images=[]

  for i,f in enumerate(files):
    p = Image.fromarray(labels[i]) 
    p.save(join(save_path, (names[i]+"_wshed"+".tif")))
    #images.append(Image.open(f))

  with ZipFile('output.zip', 'w') as zipObj:
    for folderName, subfolders, filenames in os.walk(save_path):
        for filename in filenames:
            filePath = join(folderName, filename)
            zipObj.write(filePath, os.path.basename(filePath))

  print("post processing done")
  return [labels, names, images, cdp_smooth, ndp_square, binmask, seeds]

print("Initialization done")

print(f'Pytorch version : {torch.__version__}')
print(f'torchvision version : {torchvision.__version__}')

ndp_model = load_model(UNet("ndp"),sys.argv[1]+"ndp_model.pt")
cdp_model = load_model(UNet("cdp"),sys.argv[1]+"cdp_model.pt")

prediction = predict(sys.argv[2], sys.argv[3], 0.15, 0.25)
