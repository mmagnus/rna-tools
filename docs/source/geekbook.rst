Geekbook & rna-tools
============================================================

G33KB00K - 🤓 eXtreme eXtendable note taking system for nerds/geeks (including scientists!) (= beautiful html generator of your markdown-based notes) docs: http://geekbook.rtfd.io

Marcin Magnus (mmagnus) & Pietro Boccaletto (akaped)

The code of the project can be found at GitHub (https://github.com/mmagnus/geekbook).

A neat way how to combine Emacs/Atom/Sublime editor + Markdown Syntax + Git + Html engine (bootstrap/python) to get the best notes-talking experience ever. Highly customizable with plugins written in Python. What's the most important, under the hood it's just a set of Markdown files.. you can do with them whatever you want, e.g. you can Pandoc (http://pandoc.org/epub.html) them to epub (that's origin of "book" part of the name).

Draw VARNA-based image of RNA secondary structure
------------------------------------------------------------

Type::

  <pre>[ss:rna]
  UUUCUGUAUAAUGCCGAUAAUAAGGUUCGGCAGUUUCUACCAAACAGCCGUAAACUGUUUGACUACAGUAA
  ((.(((((...((((((.........))))))........(((((((.......)))))))..))))).))
  </pre>

.. warning :: Keep exactly the same syntax as in the example above and below.

The syntax::

     <pre>
     [ss:/name of your seq/]
     /seq/
     /ss/
     </pre>
     # ^ not <pre/> nor <pre>. Keep a new line after this syntax. So don't do:
     </pre>
     <pre>

     but

     </pre>

     <pre>
     # ^ this could be fixed at some point

.. warning :: This plugin will change your Markdown file, so make sure that your editor will detect this change and ask you to reload the file!

to get a VARNA-drawn image of secondary structure.

.. image :: ../pngs/XQxIC1CrCP.gif
