c = g++ -Wall -std=c++17 -static -lm $< -o $*

%: %.cpp
	$c -g

t%: %.cpp
	$c -O2
	@for i in $*.in*; do echo "\n== $$i ==" && $(mtime) ./$* < $$i; done
