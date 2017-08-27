
===============================================================================
rasbhari 1.4.0 ReadMe by L. Hahn, August 2017, G.-A.-University of Goettingen
===============================================================================

CONTENT

#0 About rasbhari 
#1 Extract tar-gz-file
#2 Requirements
#3 Compilation
#4 Run RasBhari
#5 Authors and location
#6 Copyright
#7 Thanks!
#8 How to cite

-------------------------------------------------------------------------------
0) About rasbhari
	(Rapid Approach for Seed optimization Based on a Hill-climbing 
	Algorithm that is Repeated Iteratively)

	rasbhari is a program that creates patterns/seeds for a user chosen 
	configuration for number of patterns, their weight and their length.

	If the user wishes, rasbhari optimizes the variance or overlap
	complexity of a patternset, which were generated automatically or were 
	read from a submitted patternfile.

	It is also possible to calculate the sensitivity of a patternset and
	tries to optimize the patternset to a higher sensitivity.

-------------------------------------------------------------------------------
1) Extract zip-file

	Open the file containing folder and press right mouse-button, click
	on 'Extract here...'.

	   -or-	
	
	Use the command line and move to the containing folder. Use the 
	following command to extract the archive:
		$ tar xzfv rasbhari.tar.gz

-------------------------------------------------------------------------------
2) Requirements

	To compile and run rasbhari an C++ IDE like GCC or Visual Studio / a 
	C++ compiler and the latest revision of C++ standard(2011), also known
	as C++11 standard, is needed.

-------------------------------------------------------------------------------
3) Compilation

	To compile rasbhari on Unix/Linux systems/derivatives, use 'make' for
	the compilationafter changing into the RasBhari folder:

		$ make

	   -or-

	If you want to, or do not have make, you can manually compile rasbhari:

		$ g++ -std=c++11 -O3 -Wall main.cpp rasbopt.cpp rasbimp.cpp rasbhari.cpp
					rasbcomp.cpp sensmem.cpp speedsens.cpp patternset.cpp 
					pattern.cpp -o rasbhari

	If you do not use the GCC-compiler, the term 'g++' might change!

-------------------------------------------------------------------------------
4) Run RasBhari

	To run RasBhari please have a look at the command-help-file!
	An example program launch could be:

		./rasbhari -m 10 -w 8 -d 6-15 -H 64 --permut 25000

-------------------------------------------------------------------------------
5) Authors and location

	This software/program was created by Lars Hahn at the George-August-
	Universtiy of Goettingen.
	
	----------

	The ideas are from Burkhard Morgenstern, Lars Hahn and Chris-André
	Leimeister at the George-August-University of Goettingen.

	The following functions of the speedsens.cpp file:

		- inline long long BIN_REVERSED_TO_INT2(char *s)
		- double MULTIPLE_SENSITIVITY2(char** SEEDS, int NO_SEEDS, long long N, double P)

	after the line with the speed comment, were created by Silvana and 
	Lucian Ilie in the software called 'SpEED'.
	Therefore please have a look at: 

	Lucian Ilie, Silvana Ilie, and Anahita M. Bigvand. SpEED: fast
	computation of sensitive spaced seeds. Bioinformatics, 27:2433–2434, 2011.

-------------------------------------------------------------------------------
6) Copyright

	This software is created under the Terms of the GNU GENERAL PUBLIC LICENSE.
	Therefore please have a look at the 'Copying' file!.

-------------------------------------------------------------------------------
7) Thanks!

	A special thanks for Laurent Noé for pointing out a similarity between
	the overlap complexity from Ilie & Ilie and the variance, and also for
	some debugging for the first public run!

-------------------------------------------------------------------------------
8) How to cite
	If you are using rasbhari for your own approaces, please cite us:

	| Hahn L, Leimeister C-A, Ounit R, Lonardi S, Morgenstern B (2016)
	| rasbhari: Optimizing Spaced Seeds for Database Searching, Read Mapping and Alignment-Free Sequence Comparison
	| PLoS Comput Biol 12(10):e1005107. doi:10.1371/journal.pcbi.1005107