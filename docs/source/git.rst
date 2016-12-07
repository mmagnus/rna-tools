Git Quickref
=========================================

.. image:: ../pngs/Git-Logo-1788C.png
	   
Git is a version control system that is used for software development that helps you to keep track of versions of your program. To start using git you have to know only these two commands below. If you want to contribute to the package you need a few more, but it's not important right now :-)

To get the package for the first time on your
computer run::

  $ git clone git@github.com:mmagnus/rna-pdb-tools.git

and if you want to update the package later run::

  $ git pull # be in the folder like ~/src/rna-pdb-tools/ <here>

if you see something like this::

  $ git pull
  Already up-to-date.  

it means that your version of the package is up to date, congrats! :-)

If you see something like this::

	 $ git pull
	remote: Counting objects: 3, done.
	remote: Compressing objects: 100% (1/1), done.
	remote: Total 3 (delta 2), reused 3 (delta 2), pack-reused 0
	Unpacking objects: 100% (3/3), done.
	From github.com:mmagnus/rna-pdb-tools
	  69c4ee3..7f90739  master     -> origin/master
	Updating 69c4ee3..7f90739
	Fast-forward
	install_links_bin.sh | 1 + 
	1 file changed, 1 insertion(+)

it means that there is a small change in ``install_links_bin.sh`` and you are up to date, congrats as well! You might need to run ``./install_links_bin.sh`` to "install" new tools that were added to the packages (if this is the case). If you get any error then talk to me ``magnus@genesilico.pl``.
