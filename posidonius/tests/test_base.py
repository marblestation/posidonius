import os
import inspect
import unittest

class TestBase(unittest.TestCase):

    def setUp(self):
        self.current_filename, ignore = os.path.splitext(os.path.basename(__file__)) # Filename without extension
        self.current_dirname = os.path.dirname(inspect.getfile(inspect.currentframe()))

    def tearDown(self):
        pass
