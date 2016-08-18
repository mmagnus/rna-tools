# different sequence of atom
echo
echo 'test #1'
./diffpdb.py 'test_data/1/1kxkA.pdb' 'test_data/1/1kxkA_M1.pdb'
# one atom different
echo
echo 'test #2'
./diffpdb.py 'test_data/2/str1.pdb' 'test_data/2/str2.pdb'

echo
echo 'test #2 + htmlout'
./diffpdb.py --htmlout 'test_data/2/str1.pdb' 'test_data/2/str2.pdb' > output/2_str1_vs_str2.html
# identical
echo
echo 'test #3'
./diffpdb.py 'test_data/3.ident/str1.pdb' 'test_data/3.ident/str2.pdb'


