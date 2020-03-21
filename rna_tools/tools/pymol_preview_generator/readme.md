# PyMOL Preview Generator [OSX]

![](demo.gif)

Video with an explanation how to set it up:

<a href="https://www.youtube.com/watch?v=TdVVfPHlr8U"><img src="demo.png"></a>

1:26 Installation
4:05 Shell tool
5:52 Quick Action with Automator

Installation:

    # install brew [if you don't have it]
    /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install.sh)"

    $ brew install git # install git
    $ git clone https://github.com/mmagnus/rna-tools.git # get rna-tools [only if you don't have it]

    $ brew install fileicon
    $ brew install imagemagick
    $ brew install pymol
    
Automator:

- Workflow receives current: files
- Select Image that you want, we can use QuickLook
- Color orange
- Action: Run shell script, drag and drop this action to the Main Panel
- Copy paste the code below for this action
- Pass input `as arguments` (!)
- Save action as `PyMOL Preview`

Script for Automator:

    for f in "$@"
    do
        echo "$f"
        <PATH TO THE SCRIPT RNA-TOOLS>/rna-tools/rna_tools/tools/pymol_preview_generator/pymol_preview_generator.py $f
    done

for me:

    for f in "$@"
    do
        echo "$f"
        /Users/magnus/work-src/rna-tools/rna_tools/tools/pymol_preview_generator/pymol_preview_generator.py $f
    done

