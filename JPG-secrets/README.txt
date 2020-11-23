I made this an afternoon after experimenting with the idea of hiding a message
in an image. An image can be represented as a series of pixels, and each pixel
can be represented as three numbers between 0 and 255. We can take a message in
letter, convert it to a binary string, and then hide that binary string in the
image by modifing each pixel r, g or b value to the nearest odd or even number.
The image will not appear altered to the naked eye. Thats the theory at least,
I havent made a program to retrieve the message back from the image yet.
