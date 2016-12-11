import unittest
import subprocess

class test_sh(unittest.TestCase):
    """Super-dirty test of test.sh"""
    def test(self):
        cmd = './test.sh'
        p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        err = p.stderr.read()
        self.assertEqual(err, '')

if __name__ == '__main__':
    unittest.main()
