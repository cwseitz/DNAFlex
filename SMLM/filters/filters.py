import numpy as np
import math
from skimage.util import img_as_ubyte
from skimage.filters import gaussian
from skimage.morphology import disk, erosion, dilation
from skimage.morphology import disk
from skimage.filters import median
from skimage.morphology import disk
from skimage.filters import rank
from scipy import signal
from skimage import img_as_ubyte

def blur(image, sigma):
	image = gaussian(image, sigma)
	return image

def mean(image, disk_radius):
	selem = disk(disk_radius)
	image = rank.mean(image, selem=selem)
	return image

def median(image, disk_radius):
	selem = disk(disk_radius)
	image = median(image, selem=selem)
	return image

def minimum(image, disk_radius):
	selem = disk(disk_radius)
	image = erosion(image, selem=selem)
	image = dilation(image, selem=selem)
	return image

def boxcar(image, width):
	def normalize(array):
		array = [x/sum(array) for x in array]
		return array
	def zero(filtered_image_array, width, lnoise):

		filtered = filtered_image_array
		count = 0

		lzero = int(max(width,math.ceil(5*lnoise)))
		#TOP
		for row in range(0, lzero):
			for element in range(0, len(filtered[row])):
				filtered[row, element] = 0
		#BOTTOM
		for row in range(-lzero, -0):
			for element in range(0, len(filtered[row])):
				filtered[row, element] = 0
		#LEFT
		for row in range(0, len(filtered)):
			for element_index in range(0, lzero):
				filtered[row, element_index] = 0
		#RIGHT
		for row in range(0, len(filtered)):
			for element_index in range(-lzero, -0):
				filtered[row, element_index] = 0
		#set negative pixels in filtered to 0
		for row in range(0, len(filtered)):
			for element_index in range(0, len(filtered[row])):
				if filtered[row, element_index] < 0:
					count +=1
					filtered[row, element_index] = 0

		return filtered

	#normalize image
	image = image/255

	#build boxcar kernel
	length = len(range(-1*int(round(width)), int(round(width)) + 1))
	boxcar_kernel = [float(1) for i in range(0, length)]
	boxcar_kernel = normalize(boxcar_kernel)

	#convert 1d to 2d
	boxcar_kernel = np.reshape(boxcar_kernel, (1, len(boxcar_kernel)))
	filtered = signal.convolve2d(image.transpose(),boxcar_kernel.transpose(),'same')
	filtered = signal.convolve2d(filtered.transpose(), boxcar_kernel.transpose(),'same')
	filtered = img_as_ubyte(image-filtered)

	return filtered
