# Spotifier

## Prepare image
Open the image:
![](imgs/200618_Screen_Shot_2020-06-18_at_2.08.30_PM.png)

Image -> Adjustmenst -> Black & white

Inverse color (cmd+i) 

Adjust Levels (cmd+l), select black backgroud with the black pippette and white for the white pippette.
![](imgs/200618_Screen_Shot_2020-06-18_at_2.08.56_PM.png)
![](imgs/200618_Screen_Shot_2020-06-18_at_2.09.17_PM.png)
![](imgs/200618_Screen_Shot_2020-06-18_at_2.11.34_PM.png)

Crop image and Edit -> Image Rotation -> Flip Canvas Horizontal.

## Adjust the image to the template
Open template.psd and drag and drop plate photo
![](imgs/200618_Screen_Shot_2020-06-18_at_2.13.27_PM.png)

Lower opacity for the imported image, to aroud 30% and 
![](imgs/200618_Screen_Shot_2020-06-18_at_2.14.06_PM.png)

Use Move tool and Free transform fit the image to the template.
![](imgs/200618_Screen_Shot_2020-06-18_at_2.16.09_PM.png)

Use Move tool and Free transform to move "Rectangle" to cover the plate.
![](imgs/200618_Screen_Shot_2020-06-18_at_2.23.08_PM.png)

Switch off Backgroud layer, set Opacity to 100.
![](imgs/200618_Screen_Shot_2020-06-18_at_2.24.11_PM.png)

Save as JPG, e.g., s02_30.jpg.
![](imgs/200618_Screen_Shot_2020-06-18_at_2.25.17_PM.png)

## Prepare mapping file
Open a text editor and prepare a file used to map dots into figure. 

![](imgs/200618_Screen_Shot_2020-06-18_at_2.33.10_PM.png)

Run the program:

    python spotifier.py testdata/02/s02_30.jpg -t -m testdata/02/map.txt

The results should be like this:

![](imgs/200618_Screen_Shot_2020-06-18_at_3.03.08_PM.png)

and the file `s02_30_spots.png` should be created in the folder next to the input file (in this case `testdata/s02/s02_30_spots.png`)


