One of the assigments in High performance computing I was given was an image
denoising (or smoothing) program. The task was to first make the program in
serial C for a black and white image, and then make the same in parallell
code on a color image. As you may see from the results, I failed the second
task. However, the glitches in the processed image showcase how the image was
divided perfectly between the cpu cores, so that the number of edge pixels were
minimized. The division would in theory be optimized for memory bandwidth.

Note that the program will not properly compile and execute without the
appropriate libjpeg dll and libraries.
