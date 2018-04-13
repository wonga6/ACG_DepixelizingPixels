# ACG_DepixelizingPixels

Final Project for Advaced Computer Graphics at Rensselaer Polytechnic Institute

C++ Implementation of the paper "Depixelizing Pixel Art" by Johannes Kopf and Dani Lischinski:
https://johanneskopf.de/publications/pixelart/paper/pixel.pdf
Notes:
- not a complete implementation of the paper - last section worked on in paper was 3.3
    - shapes aren't combined or smoothed out
- goal (accomplished) was to create an svg image from a pixel image that could be edited in a vector drawing program

To run:

g++ -g main.cpp image.cpp edge.cpp node.cpp -o test.out -Wall -std=c++11

./test.out ../TestImages/Image.ppm


References:
[1] 2015. Color Conversion. Equasys. Web, http://www.equasys.de/colorconversion.html  
[2] Cutler, B/ 2015. Homework 9. Cs.rpi.edu. Web.  
[3] Kopf, J. & Lischinski, D., 2011. Depixeling Pixel Art. ACM Transactions on Graphics (Procedings of SIGGRAPH 2011),30, 4, 99:1 - 99:8.  
[4] Selinger, P., 2003. Potrace: a polygon-based tracing algorithm. http://potrace.sourceforge.net  
[5] Turney, M., 2010. simple-svg. Github repository, https://github.com/adishavit/simple-svg  
[6] Vemula, A. & Yeddu, V., 2014. Pixel Art. Github repository, https://github.com/vvanirudh/Pixel-Art  
[7] Wala, R.J., 2016. depixelize. Github repository, https://github.com/rjalfa/depixelize
