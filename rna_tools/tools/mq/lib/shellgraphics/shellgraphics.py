#!/usr/bin/python
"""ShellGraphics - a library of some simple graphics tricks in terminal

Run `ShellGraphics.py` or open `demo.txt` (the file was created with `python ShellGraphics.py > demo.txt`)

 ____  _          _ _  ____                 _     _          
/ ___|| |__   ___| | |/ ___|_ __ __ _ _ __ | |__ (_) ___ ___ 
\___ \| '_ \ / _ \ | | |  _| '__/ _` | '_ \| '_ \| |/ __/ __|
 ___) | | | |  __/ | | |_| | | | (_| | |_) | | | | | (__\__ \
|____/|_| |_|\___|_|_|\____|_|  \__,_| .__/|_| |_|_|\___|___/
                                     |_|                     

 mqaprna.py

####################################################################################################
# mqaprna.py                                                                                       #
####################################################################################################
****************************************************************************************************
* mqaprna.py                                                                                       *
****************************************************************************************************

 mqaprna.py

 Ignore pdb fn        off
 Seq ss fn            off
 Native pdb           off
 Output csv           test_m0.csv
 Directory            test_data/2/
 Multiprocessing      off
 AnalyzeGeometry      off
 SSAgreement          off
 to ignore []
 files to analyze: 3
----------------------------------------------------------------------------------------------------
 methods: ['RASP', 'SimRNA', 'ClashScore', 'FARNA', 'NAST_pyro']
----------------------------------------------------------------------------------------------------

etc.
"""
import sys
import time
import subprocess

terminal_length = 80
color_mode = True

class Title(object):
    def __init__(self):
        print("#" * terminal_length)
        title = "mqaprna.py v0.1"
        print("#", title, " " * (terminal_length  - 5 - len(title)), '#')
        print("#" * terminal_length)

class Table(object):
    column_length = 30
    def __init__(self, headers):
        self.headers = headers

    def show_header(self):
        print(self.headers[0].ljust(self.column_length), end=' ')
        for h in self.headers[1:]:
            print(h.rjust(self.column_length), end=' ')
        print()
        print('-' * terminal_length)

    def show_row(self, values):
        print(values[0].ljust(self.column_length), end=' ')
        for h in values[1:]:
            print(h.rjust(self.column_length), end=' ')
        print()

    def save_csv(self):
        pass


def table_demo():
    t = Table(['fn','SimRNA', 'FARNA', 'X'])
    t.show_header()
    t.show_row(['1xjrA_output3-000142_AA.pdb', '-35.861', '44.043', '111'])
    t.show_row(['1xjrA_output3-000142_AA.pdb', '-34.861', '322.12', '443'])
    phr()


def phr(char='-'):
    """Print ------- where '-' is set by char"""
    print(char * terminal_length)


def phr_text(text, verbose=True):
    """Print hr if CHAR is '*' you will get:
    ************************************************************** <text> *"""
    len_of_hr = terminal_length
    CHAR =  '-'
    if verbose: 
        text=str(text)
        st=CHAR * (len_of_hr - 3 - len(text)) + ' ' + text + ' ' + CHAR # 3 = ' 'test' '* = 2x' ' + 1x*
        print(st)
        return True
    else:
        return False


def pbanner_simply(text):
    """
    """
    pbr()
    if color_mode:
        print('\033[1;34m' + text + '\033[1;m')
    else:
        print(text)
    pbr()


def pbanner(text, char="#", verbose=True):
    """Print

    - text, like: SUPER PROGRAM ver.01
    - len, by default it equals to LENGTH_OF_THE_SCREEN
    - verbose
    ---------

    DO:
    - print banner, like:
    ##################################################
    # SUPER PROGRAM ver.01                           #
    ##################################################

    or use '*'
    ---------
    RETURN:
    - True if ok
    """
    len = terminal_length
    if verbose:
        print(char * len + '\n' + char +' '+ text.ljust(len-3)+ char  +'\n'+ char * len)
    return True


def pprogress_line(step, amount,text,length_the_bar=terminal_length, verbose=True):
    """Generate a progress bar::

      [===================================================X----] # 36412 ; 94.0878552972 % ; 36500

    :params:

     - step, c(ounter) in loop,
     - amount, to count a percent,
     - length_the_bar of whole progress bar (use to scale)

    how to use it::

        c=0
        len_of_iter=len(list)
        for i in list:
            c+=1
            toolboX.progress_line(c,len_of_iter,verbose=True)
    
    :return: string
    """
    percent=float(step)/amount * 100
    #
    #  length_the_bar 
    #
    #length_the_bar=length_the_bar#-15
    #
    done=int((float(percent)*length_the_bar)/100) # 1.4 % -> 1
    #
    #  2% calosci - 100%
    #  x (ile kresek) - 50 (gdy dlugosc bara wynosi jeno 50); czyli wyjdzie o polowe mnieij
    #  dla procent 2% tylko jeden step done '[=X ...'
    #
    #for done in range(1,100): to test
    #print '#',step,';',percent, '%'            
    if verbose: print('['+'-'*done+''+' '*(length_the_bar-done)+']',' ',str(step).rjust(4), str(round(percent,2)).ljust(4), '%',amount, text)
    # return True


def pformat_line(name, value='',char="#", verbose=True, wide=terminal_length):
    """
    GET:
    - name,
    - value,
    - verbose,
    - wide, length of field with name, len from the beginning of name to ':'


    DO:
    print # STEP              : 20

    pformat_line('magnus','ok')
    # Magnus                                                            o
    """

    value=str(value)

    if verbose: 
        if value:
            #print "# "+name.capitalize().ljust(wide)+' : ',value # #EXISTS               : Yes
            if type(value)==type([1]):
                print(char + ' ' + name.capitalize() +' ', value)
            else:
                print(char + ' ' + name.capitalize()+' '* (LENGTH_OF_THE_SCREEN - len(value) - len(name)-2) +value) #EXISTS               : Yes
        else:
            print(char + ' ' + name.capitalize().ljust(wide))               #EXISTS              


def pin_red(text, newline = True):
    """
    http://www.siafoo.net/snippet/88
    """
    if color_mode:
        if newline:
            print('\033[1;31m' + text + '\033[1;m')
        else:
            print('\033[1;31m' + text + '\033[1;m', end=' ')        
    else:
        print(text)
    
def pin_green(text, newline = True):
    """
    http://www.siafoo.net/snippet/88
    """
    if color_mode:
        if newline:
            print('\033[1;32m' + text + '\033[1;m')
        else:
            print('\033[1;32m' + text + '\033[1;m', end=' ')
    else:
        print(text)

def pin_blue(text, newline = True):
    """
    http://www.siafoo.net/snippet/88
    """
    if color_mode:
        if newline:
            print('\033[1;34m' + text + '\033[1;m')
        else:
            print('\033[1;34m' + text + '\033[1;m', end=' ')
    else:
        print(text)


def pin_red_and_blue(text, text2, newline = True):
    """

    inspired: http://www.siafoo.net/snippet/88
    """
    if color_mode:
        if newline:
            print('\033[1;31m' + text + '\033[1;m' + '\033[1;34m' + text2 + '\033[1;m')
        else:
            print('\033[1;31m' + text + '\033[1;m' + '\033[1;34m' + text2 + '\033[1;m')
    else:
        print(text + ' ' + text)


def ph_center(text, char='-'):
    if len(text) % 2:
        half = (terminal_length - len(text) - 1) / 2
        print(char * (half - 1), text, char * half)
    else:
        half = (terminal_length - len(text) - 2) / 2
        print(char * half, text, char * half)


def ph_center2(text):
    """
    ================================================================
                            Output                              
    ----------------------------------------------------------------
    """
    phr('=')
    ph_center(text, char=' ')
    phr('-')


def ptitle(text):
        """Show a huge title using `figlet`.

        .. warning :: `figlet` has to be installed in system
        """
        subprocess.check_call(['figlet', text])


def poption(name, value):
    """Print option: name and value. If value is True or 'on' or any string/number then print 'on' in green otherwise print 'off' in red.
    """
    if (value == 'off' or value == False or not value) and color_mode:
        print('', name.ljust(20), '\033[1;31moff\033[1;m')  # red
    elif (value == 'on' or value == True) and color_mode:
        print('', name.ljust(20), '\033[1;32mon\033[1;m')  # green
    elif value and color_mode:
        print('', name.ljust(20), '\033[1;32m' + value + '\033[1;m')  # green
    else:
        print('', name.ljust(20), value)


def poptions(dictionary):
    """Get dict and print as options
    """
    print('--------------------------  ----------------------------------------------------')
    for d in list(dictionary.keys()):
        print(d.ljust(27), end=' ') 
        value = dictionary[d]        
        if (value == 'off' or value == False or not value) and color_mode:
            print('\033[1;31moff\033[1;m')  # red
        elif (value == 'on' or value == True) and color_mode:
            print('\033[1;32mon\033[1;m')  # green
        elif value and color_mode:
            print('\033[1;32m' + value + '\033[1;m')  # green
        else:
            if value == None:
                print('off')
            elif value == False:
                print('off')
            elif value == True:
                print('on')
            else:
                print(value)


def pbr():
        print()


def psub(text):
    """Print ' <text>'
    """
    print('', str(text))

def p(text):
    """Just print '<text>'
    """
    print(str(text))


def start(): pass

if __name__ == '__main__':
        ptitle('ShellGraphics')
        pbanner_simply('mqaprna.py')
        pbanner('mqaprna.py')
        pbanner('mqaprna.py', '*')
        phr_text('SimRNA')

        for i in range(1,6):
                pprogress_line(i,5, 'text')
                sys.stdout.write("\033[F") # Cursor up one line - works in oneline!
                time.sleep(1)

        Title()
        table_demo()
        pformat_line('foo foo')
        pin_red('foo')
        pin_blue('foo')
        pin_red_and_blue('foo','foo2')
        ph_center('Options')
        p('help')
        psub('logout - to log out')
        ph_center('Options:')
        ph_center2('Option')
        pbr()
        poption('native pdb', 'on') 
        poption('native pdb', 'off') 
        pbr()
        poption('native pdb', True) 
        poption('native pdb', False) 
        pbr()
        poption('native pdb', 'any text')
        poption('native pdb', bool('any text'))
        poption('native pdb', '') 
