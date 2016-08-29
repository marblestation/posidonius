#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Script that is a Makefile-like, but in python, and for fortran90

import glob
import sys
import os
import string
import subprocess

LOG_NAME = 'compilation.log'

#COMPILATOR = "gfortran"

#DEBUG_OPTIONS = "-pedantic-errors -Wall -Wconversion -Wextra -Wunreachable-code -fbacktrace" + \
  #" -g3 -fbounds-check -O0" + \
  #" -fstack-protector-all -fno-automatic -Wuninitialized -ftrapv -fno-automatic -fimplicit-none"
## commented options : -ffpe-trap=invalid,zero,overflow,underflow because most of the time, this is not a bug.
#OPTIMIZATIONS = "-O3 -ffast-math -pipe -finit-real=nan"
#TEST_OPTIONS = "-O3 -pipe -finit-real=nan"
#GDB_OPTIONS = "-g3"
#PROFILING_OPTIONS = "-g -pg"

COMPILATOR = "ifort"

DEBUG_OPTIONS = "-g -check all -fpe0 -warn -traceback -debug extended"
OPTIMIZATIONS = "-O3 -Bstatic"
TEST_OPTIONS = OPTIMIZATIONS
GDB_OPTIONS = "-g"
PROFILING_OPTIONS = "-g -pg"

class sourceFile(object):
  """Define an object linked to a fortran 90 source code that will
  store dependencies, including modules included, used or included

  Parameters :
  filename : the name of an existing source code filename (for example "main.f90")
  name='default' : the name we want for a program (if it is a program, if not, there will be no use)
  isProgram=False : is this source file a main source for a binary we want, or just a module or sub-program?

  Attribute :
  self.filename : the filename of the source code file.
  self.defined : a list of procedures defined in the fortran source code
  self.used : a list of modules that are used by the code
  self.included : a list of things included in the code
  self.isCompiled : a boolean to say if the source file has already been compiled.
  self.dependencies : a list of object (*.o) filenames we need to compile
  self.isProgram : a boolean to say if we want to have a binary, or if it's just a module or a subprogram
  self.name : the name we want for the binary file if it is a program. This name is by default the filename without the extension

  Methods :
  .compile() : Compile the current program and all the required dependencies

  """

  ## We define a dictionary where we store, for each name of a module,
  # the link toward the object file where it is defined
  findModule = {}

  # We define a dictionary to make the correspondance between a source filename and the object associated
  findSource = {}
  #~ COMPILATOR = "ifort"
  #~ OPTIONS = "-vec-report0 -i-dynamic -mcmodel=medium -shared-intel -L/usr/lib64/atlas -llapack"

  COMPILATOR = "gfortran"
  OPTIONS = "-O3 -march=native"

  def __init__(self, filename, name=None, isProgram=False, extra_files=None):
    """Will check everything that is included in the source code
    and initialize the object"""
    self.filename = filename

    if (name == None):
      self.name = self.filename.rstrip(".f90")
    else:
      self.name = str(name)

    # For extra modules that can't be retrieved manually. This can be C modules or weird files
    # Adding manually the modules that fuck everything up, namely ODEPACK !
    if (extra_files != None):
      # If the source file is newer than the object file, we need to compile it
      object_files = ["%s.o" % os.path.splitext(filename)[0] for filename in extra_files]
      self.extra = " ".join(object_files)

      for (source_file, object_file) in zip(extra_files, object_files):
        (name, ext) = os.path.splitext(source_file)

        if (os.path.isfile(object_file)):
          # If source file is newer, we compile
          isCompilation = os.path.getmtime(source_file) > os.path.getmtime(object_file)

        else:
          # If the object file do not exist
          isCompilation = True

        if isCompilation:

          options = sourceFile.OPTIONS

          if (ext == '.f90'):
            commande = sourceFile.COMPILATOR+" "+sourceFile.OPTIONS+" -c "+source_file
          elif (ext == '.c'):
            commande = "gcc -c "+source_file
          else:
            raise ValueError('Unkown extension for source file: %s' % ext)

          process = subprocess.Popen(commande, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)

          (process_stdout, process_stderr) = process.communicate()

          print("Compiling "+source_file+"...")
          returnCode = process.poll()


          # if returnCode is not 0, then there was a problem
          if (returnCode != 0):
            # We write compilation errors in the following file.
            f = open(LOG_NAME,'w')
            f.write(process_stderr)
            f.close()

            print("Compilation error, see '%s'" % LOG_NAME)
            LogPostProcessing()
            sys.exit(1)
          else:
            if (len(process_stderr) != 0):
              # We write compilation errors in the following file.
              f = open(LOG_NAME,'a')
              f.write(process_stderr)
              f.close()

    else:
      self.extra = ""

    # By default, nothing is a program
    self.isProgram = isProgram


    self.isCompiled = False

    if (not(self.isProgram)):
      # If the source file is newer than the object file, we need to compile it
      object_file = "%s.o" % self.name
      if (os.path.isfile(object_file) and os.path.isfile("%s.mod" % self.name)):
        self.toBeCompiled = os.path.getmtime(self.filename) > os.path.getmtime(object_file)
      else:
        # If the object file do not exist
        self.toBeCompiled = True
    else:
      self.toBeCompiled = True # We always compile the programs

    (self.defined, self.used, self.included) = self.__getModules()

    for module in self.defined:
      sourceFile.findModule[module] = self

    sourceFile.findSource[self.filename] = self

  @classmethod
  def setCompilingOptions(cls, options):
    """method that set the 'OPTIONS' value.

    Parameter:
    options : a string with the compilation options for the binary construction
    """

    cls.OPTIONS = options

  @classmethod
  def setCompilator(cls, compilator):
    """method that set the 'COMPILATOR' value.

    Parameter:
    options : a string with the compilator you want to use
    """

    cls.COMPILATOR = compilator

  def __getModules(self):
    """returns a tuple containing the list of defined modules and
    the list of used modules of a fortran source file

    Return
    defined : a list of procedures defined in the fortran source code
    used : a list of modules that are used by the code
    included : a list of things included in the code
    """

    f=open(self.filename,'r')
    lines = f.readlines()
    f.close()

    defined=[]
    used=[]
    included=[]

    for lsave in lines:
      l=string.expandtabs(string.lower(lsave)[:-1],1)
      words=string.split(string.lstrip(l))
      if len(words) > 0:
        if words[0] == 'use':
          used.append(string.split(words[1],',')[0])
        if words[0] == 'module':
          if len (words) == 2 or words[1] != "procedure":
            defined.append(words[1])
        if words[0] == 'include':
          newstring = string.replace(words[1],'\'','')
          newstring = string.replace(newstring,'\"','')
          included.append(newstring)
      l=string.expandtabs(lsave[:-1],1)
      words=string.split(string.lstrip(l))
      if len(words) > 0:
        if words[0] == '#include':
          newstring = string.replace(words[1],'\"','')
          included.append(newstring)

# We delete all dependencies that are present several number of times.
    used = list(set(used))

    return defined,used,included

  def __getFirstOrderDependence(self):
    """return a list of *.o files we need to compiled. They are
    extracted from direct used files, that why we call them
    "first order" dependance."""

    dependances = []
    for mod in self.used:
      try:
        source = sourceFile.findModule[mod]
      except:
        print("Error: Unable to locate the module '"+mod+"'")
      obj = string.replace(source.filename,'.f90','.o')
      dependances.append(obj)

    return dependances

  def setProgram(self, boolean):
    """method to defined the current object as a program, that is, if
    we want to have a binary from this source"""

    self.isProgram = boolean

  def compile(self, parent_dependencies=[]):
    """method that check dependencies and try to compile

    Parameter :
    parent_dependencies=[] : list that store all the parent dependencies of the current file, namely, all the module
    that MUST NOT be used inside the current module, else we will have an infinite loop.

    Beware, there MUST NOT be inter dependant modules.
    """

    parent_dependencies.append(self.name)

    if not(self.isCompiled):
      # We only store, for the moment, the first order dependencies.
      self.dependencies = self.__getFirstOrderDependence()

      # We store links towards all the object sourceFile that defined the modules we are interested in.
      module_sources = []
      for module in self.used:
        module_sources.append(sourceFile.findModule[module])

      # For each object, we check if there is loop call of modules. If that's the case, return an error
      for name in self.used:
        if name in parent_dependencies:
          error_message = "The module '"+self.name+"' try to use the module '"+name+"' that already "+ \
                          "use the module '"+self.name+"'. So there is an infinite loop that is not correct."
          raise NameError(error_message)

      # For each object, we compile it if it's not already the case.
      for source in module_sources:
        # We must launch source.compile() even if the file will not be compilated, just to get the right dependencies.
        if not(source.isCompiled):
          # the list() is here to ensure not to have a pointer and share the list. If not, the list of parent_dependencies will
          # not be correct and contains all the previous parent dependencies.
          source.compile(list(parent_dependencies))

          # We compile the current file if any of the dependencies has been effectively compiled
          if source.toBeCompiled:
            self.toBeCompiled = True

      # We complete the dependencies list now that all used modules
      # have been compiled, they must have a complete list of
      # their own dependencies.
      for source in module_sources:
        self.dependencies.extend(source.dependencies)

      # We delete all dependencies that are present several number of times.
      self.dependencies = list(set(self.dependencies))

      if self.toBeCompiled:
        # Now that all the dependencies have been compiled, we compile
        # the current source file.

        if not(self.isProgram):
          commande = sourceFile.COMPILATOR+" "+sourceFile.OPTIONS+" -c "+self.filename
        else:
          commande = sourceFile.COMPILATOR+" "+sourceFile.OPTIONS+" -o "+self.name+" "+self.filename+" "+self.extra+" "+" ".join(self.dependencies)
          print(commande)

        process = subprocess.Popen(commande, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)

        (process_stdout, process_stderr) = process.communicate()

        print("Compiling "+self.filename+"...")
        returnCode = process.poll()


        # if returnCode is not 0, then there was a problem
        if (returnCode != 0):
          # We write compilation errors in the following file.
          f = open(LOG_NAME,'w')
          f.write(process_stderr)
          f.close()

          print("Compilation error, see '%s'" % LOG_NAME)
          LogPostProcessing()
          sys.exit(1)
        else:
          if (len(process_stderr) != 0):
            # We write compilation errors in the following file.
            f = open(LOG_NAME,'a')
            f.write(process_stderr)
            f.close()


      self.isCompiled = True

  def __str__(self):
    """overload the str method. As a consequence, you can print the object via print name_instance

    return : a string that represent the properties of the object
    """

    texte = "filename: "+self.filename+"\n"

    texte += "defined: "+str(self.defined)+"\n"
    texte += "used: "+str(self.used)+"\n"
    texte += "included: "+str(self.included)+"\n"

    return texte

def prepare_compilation():
  """function that will prepare the compilation by defining options and listing
  all the source files to create the corresponding objects.

  """

  sources_filename = glob.glob("*.f90")

  # We define the objects for each source file.
  sources = []
  for filename in sources_filename:
    source = sourceFile(filename)
    sources.append(source)


def compile_source(filename, name=None, extra=None):
  """
  """

  source = sourceFile(filename, name=name, isProgram=True, extra_files=extra)
  source.compile()

def run_command(commande):
  """lance une commande qui sera typiquement soit une liste, soit une
  commande seule. La fonction renvoit un tuple avec la sortie,
  l'erreur et le code de retour"""
  if (type(commande)==list):
    process = subprocess.Popen(commande, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
  elif (type(commande)==str):
    process = subprocess.Popen(commande, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
  else:
    raise TypeError("The command is neither a string nor a list.")
  (process_stdout, process_stderr) = process.communicate()
  returncode = process.poll()
  # there is .poll() or .wait() but I don't remember the difference. For some kind of things, one of the two was not working
  return (process_stdout, process_stderr, returncode)

def is_non_committed_modifs():
  """function that return as a boolean if theere is non committed modifications in the repository"""
  (stdout, stderr, returnCode) = run_command("git diff|wc -l")

  if (returnCode != 0):
    return None

  nbLines = int(stdout)

  return (nbLines != 0)

def LogPostProcessing():
  """Function to modify LOG_NAME file in various conditions.
  Especially, supress warnings due to opkd package that we cannot modify
  since this is an unfortunate black box."""

  if ignoreOpkdWarnings:
    if os.path.isfile(LOG_NAME):
      objectFile = open(LOG_NAME, 'r')
      lines = objectFile.readlines()
      objectFile.close()

      new_lines = []
      i = 0
      while (i<len(lines)):
        line = lines[i]

        if (line.startswith("opkd")):
          for j in range(5):
            del(lines[i])
        else:
          i += 1
          new_lines.append(line)

      objectFile = open(LOG_NAME, 'w')
      for line in new_lines:
        objectFile.write(line)
      objectFile.close()

def clean(exts):
  """supprime les fichiers correspondant à l'expression donnée. La fonction renvoit la sortie si ça c'est bien
  déroulée, sinon ne renvoit rien."""
  for ext in exts:

    commande = "rm *."+ext

    process = subprocess.Popen(commande, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)

    (process_stdout, process_stderr) = process.communicate()
    returnCode = process.poll()

  # If returnCode is not 0, then there was a problem (in fact, will only do something for the last one
  if (returnCode==0):
    return process_stdout
  else:
    return returnCode



#############################################
# Beginning of the code
#############################################


# Parameters
debug = False
isTest = False
isManual = False # If true, we specifie a source file name to be compiled
isMercury = True
isElement = True
isClose = True
gdb = False
profiling = False
force = False # To force the compilation of every module
ignoreOpkdWarnings = True

isProblem = False
problem_message = """AIM: Compilation of mercury programs.
By default, all of them, but one can specify one specific code to compile
and avoid the others (mercury, element, close). The compilation of dependances is automatic.
The compilation options are packed into 3 meta-options : test, debug and gdb.
The modules that haven't changed since last compilation are not compiled again.
If you want to force compilation of all modules, use the "force" option.

The script can take various arguments:
(no spaces between the key and the values, only separated by '=')
 * help : display a little help message on HOW to use various options
 * force : To force the compilation of every module even those not modified
 * name=source.f90 : To compile a specific code
 * mercury : To compile mercury only
 * element : To compile binary for outputs (element) only
 * close : To compile binary for close encounters (close) only
 * test : [%s] activate test options. Theses options are to be used when we want
  to compare the original and actual version of the code, by
  launching tests_mercury.py
 * debug : [%s] activate debug options
 * gdb : [%s] activate options for gdb
 * profiling : [%s] activate options for profiling

 Example :
 Makefile.py gdb""" % (isTest, debug, gdb, profiling)

value_message = "/!\ Warning: %s does not need any value, but you defined '%s=%s' ; value ignored."

# We get arguments from the script
for arg in sys.argv[1:]:
  try:
    (key, value) = arg.split("=")
  except:
    key = arg
    value = None
  if (key == 'debug'):
    debug = True
    if (value != None):
      print(value_message % (key, key, value))
  elif (key == 'name'):
    isManual = True
    isMercury = False
    isElement = False
    isClose = False
    source_name = value
  elif (key == 'force'):
    force = True
    if (value != None):
      print(value_message % (key, key, value))
  elif (key == 'test'):
    isTest = True
    if (value != None):
      print(value_message % (key, key, value))
  elif (key == 'mercury'):
    isMercury = True
    isElement = False
    isClose = False
    if (value != None):
      print(value_message % (key, key, value))
  elif (key == 'element'):
    isMercury = False
    isElement = True
    isClose = False
    if (value != None):
      print(value_message % (key, key, value))
  elif (key == 'close'):
    isMercury = False
    isElement = False
    isClose = True
    if (value != None):
      print(value_message % (key, key, value))
  elif (key == 'gdb'):
    gdb = True
    if (value != None):
      print(value_message % (key, key, value))
  elif (key == 'profiling'):
    profiling = True
    if (value != None):
      print(value_message % (key, key, value))
  elif (key == 'help'):
    isProblem = True
    if (value != None):
      print(value_message % (key, key, value))
  else:
    print("the key '%s' does not match" % key)
    isProblem = True

if isProblem:
  print(problem_message)
  exit()

isModifs = is_non_committed_modifs()

# We clean undesirable files. Indeed, we will compile everything everytime.
if force:
  clean(["o", "mod"])



# Before compiling, we delete the previous compilation log. Indeed, we need to append the several warnings in the same file
# But we do not want to have infos of the previous compilation in it.
if os.path.isfile(LOG_NAME):
  os.remove(LOG_NAME)

# We create the binaries
if debug:
  OPTIONS = DEBUG_OPTIONS
elif isTest:
  OPTIONS = TEST_OPTIONS
else:
  OPTIONS = OPTIMIZATIONS

if gdb:
  OPTIONS = GDB_OPTIONS

if profiling:
  OPTIONS = PROFILING_OPTIONS

sourceFile.setCompilator(COMPILATOR)
sourceFile.setCompilingOptions(OPTIONS)

prepare_compilation()

if (isManual):
  compile_source(filename=source_name)

if (isMercury):
  compile_source(filename="mercury.f90")

if (isElement):
  compile_source(filename="element.f90")

if (isClose):
  compile_source(filename="close.f90")

if (isModifs):
  print("Warning: There is non committed modifs!")

LogPostProcessing()

if os.path.isfile(LOG_NAME):
  if 'Warning' in open(LOG_NAME).read():
    print("Warnings: see '%s'" % LOG_NAME)
