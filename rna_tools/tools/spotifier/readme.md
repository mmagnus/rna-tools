# Spotifier

This script to covert pre-processed plates image into figures.
![](imgs/yWwr8DqnXT.gif)

## Options

    $ spotifier.py -h
    usage: spotifier.py [-h] [-v] [-d] [--dont-align] [-x X] [-y Y]
                        [--trim-rms TRIM_RMS] [--size SIZE] [-a]
                        map file [file ...]

    install to work on psd files: psd-tools3

    positional arguments:
      map                  map
      file                 pre-processed image(s)

    optional arguments:
      -h, --help           show this help message and exit
      -v, --verbose        be verbose
      -d, --debug          be even more verbose
      --dont-align         don't align dots
      -x X
      -y Y
      --trim-rms TRIM_RMS
      --size SIZE
      -a, --dont-annotate

![](imgs/jyvsjrgxpe.gif)

Fig. `-t`,  trim background to get nicely formatted dots.

## Prepare image
Open the image:
![](imgs/200618_Screen_Shot_2020-06-18_at_2.08.30_PM.png)

Image -> Adjustments -> Black & white

Inverse color (cmd+i) 

Adjust Levels (cmd+l), select black backgroud with the black pippette and white for the white pippette.
![](imgs/200618_Screen_Shot_2020-06-18_at_2.08.56_PM.png)
![](imgs/200618_Screen_Shot_2020-06-18_at_2.09.17_PM.png)
![](imgs/200618_Screen_Shot_2020-06-18_at_2.11.34_PM.png)

Crop image and Edit -> Image Rotation -> Flip Canvas Horizontal.

## Adjust the image to the template
Open template.psd and drag and drop plate photo
![](imgs/200618_Screen_Shot_2020-06-18_at_2.13.27_PM.png)

Lower opacity for the imported image, to around 30%, to see the template in the background.
![](imgs/200618_Screen_Shot_2020-06-18_at_2.14.06_PM.png)

Use Move tool and Free transform fit the image to the template.
![](imgs/200618_Screen_Shot_2020-06-18_at_2.16.09_PM.png)

Use Move tool and Free transform to move "Rectangle" to cover the plate.
![](imgs/200618_Screen_Shot_2020-06-18_at_2.23.08_PM.png)

Switch off the Backgroud layer, set Opacity to 100.
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

-------------------------------------------------------------------------------

TIP: If you want to introduce one extra, empty line, just add an empty line to a map file.
![](imgs/map_empty.png)

# Customization

If you want to move single dots, use Preview and just move them around, save it (if you open a JPG, you will be asked if you want to convert the file to PNG, yeah, do it, remember only to change the file name in the command, `18_X.png`).

![](imgs/fix.png)

and re-run:

	python spotifier.py testdata/02/18_X.png -t -m testdata/02/map.txt

-------------------------------------------------------------------------------

Use play with Opacity to see Backgound with the template to adjust positions of your dosts.

![](imgs/opacity.png)

# Mappings
Some ideas for your plates ;-):

    5, 6, 7
    8, 9, 10

    12, 13, 14
    15, 16, 17

    20, 21, 22
    23, 24, 25

    28, 29, 30
    31, 32, 33

    36, 37, 38
    39, 40, 41

    44, 45, 46
    47, 48, 49
    
-------------------------------------------------------------------------------

    1, 2, 3
    5, 6, 7
    8, 9, 10
    12, 13, 14
    15, 16, 17
    21, 22, 23
    24, 25, 26
    28, 29, 30
    31, 32, 33
    36, 37, 38
    39, 40, 41
    44, 45, 46
    47, 48, 49
    51, 52, 53
    53, 54, 55
    57, 58, 59
    

# Install

    pip install rna-tools
    
or 

    pip install -e git+http://github.com/mmagnus/rna-tools.git#egg=rna-tools

and install requirements:

    pip install pillow
