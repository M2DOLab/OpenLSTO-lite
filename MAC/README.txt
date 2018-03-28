# INSTALLATION FOR MAC USERS 
> (Drafted by Carolina Jauregui)

Follow these steps to properly compile and run the main.cpp file

1. Install homebrew (skip this step if you already have homebrew)

   Copy and paste the following command into the terminal prompt
   /usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"


2. Install gcc (skip this step if gcc is already installed)
   Note: this step will take a very long time (up to 3 hours)

   Copy and paste the following command into the terminal prompt
   brew install gcc


4. Check the path of gcc

   The previous command will give output that shows the path of where gcc has been installed. Look for the prefix
   flag: --prefix. You will see something like this: --prefix=/usr/local/Cellar/gcc/7.2.0
   In this case, the path is /usr/local/Cellar/gcc/7.2.0 is the path

   If you can't find the flag in the output, do the following

   Copy and paste the following command into the terminal prompt 
   brew info gcc

   The above command will show output which contains the location of where gcc has been installed. It should be
   located under the line https://gcc.gnu.org/ and look something like this: /usr/local/Cellar/gcc/7.2.0


5. Export gcc path

   Replace GCC_PATH in the command below with the path you got in the previous step and VERSION with the first 
   number of the version of gcc that was installed. 

   Copy and paste the following command into the terminal prompt
   export CXX=GCC_PATH/bin/g++-VERSION

   EXAMPLE:
   GCC_PATH = /usr/local/Cellar/gcc/7.2.0
   VERSION = 7
   export CXX=/usr/local/Cellar/gcc/7.2.0/bin/g++-7

6. Assign CXX variable in makefile

   The first line in the makefile is: CXX = 
   Copy and paste the file directory from the previous step into the first line.

   EXAMPLE:
   CXX = /usr/local/Cellar/gcc/7.2.0/bin/g++-7

7. Make and run executable

   Copy and paste the following command into the terminal prompt
   make && ./bin/a.out
