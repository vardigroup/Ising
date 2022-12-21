# UAI ISING SPECIFICATION
The following file format is used to specify Ising models. It is based on the UAI format (https://web.archive.org/web/20170425162729/http://graphmod.ics.uci.edu/uai08/FileFormat).

## Preamble
All Ising model files start with a preamble. The preamble contains metadata about the file as well as the scopes
of each interaction. The values that interactions take are given later in the function table.

### Metadata
The metadata section of the preamble starts with "ISING" to signal that the file describes an Ising model.
Then, the number of lattice sites in the model is given. On the next line, the size of the scope of each variable
is specified - for an Ising model this should be a sequence of 2's, but for a Potts model this may vary.\
Finally, the file gives the number of field interactions (the number of lattice sites i for which h<sub>i</sub> ≠ 0) 
followed by the number of pairs of lattice sites with non-zero interaction (the number of lattice site pairs (i,j) with J<sub>i,j</sub> ≠ 0). The inverse temperature (β) and field orientation (μ) may be optionally specified on the same line as the numbers of interactions. Note that if field orientation is specified, inverse temperature must be specified as well.

For example, the following is a valid metadata section of a system with 4 lattice sites and an inverse temperature of 5:
<pre>
# This is a comment
ISING
4
2 2 2 2 # This is also a comment
4 3 5
</pre>
### Interaction Scopes

Next, the file specifies the scope of each interaction. For each field interaction, the file has a line of the form "1 i" where i denotes the ith lattice site. This means that there is a field interaction affecting the ith node. It is recommended that these lines be ordered by "i" for legibility. For each pairwise interaction between lattice sites, the file has a line of the form "1 i j" where i and j are the interacting lattice sites. It is recommended that i be less than j so that the generated matrix is upper triangular. It is also recommended that these lines be ordered first by i and then by j

For example, the interaction scope section of our file might look like this:
<pre>
1 0
1 1
1 2 # I am a comment as well
1 3
2 0 1
2 0 3
2 1 2
</pre>

Between the preamble and the function table, there is an empty line.
The function table gives the values taken by the interaction functions acting on the lattice sites specified in the interaction scope section of the file. Note that the functions in the function table should be given in the same order as the functions' scopes are listed in the interaction scope section.

Each field interaction should be written as follows:
<pre>
2
 h<sub>i</sub>(-1) h<sub>i</sub>(1)
</pre>

Each pairwise interaction between two lattice sites should be written as follows:
<pre>
4
 J<sub>i,j</sub>(-1,-1) J<sub>i,j</sub>(-1,1)
 J<sub>i,j</sub>(1,-1) J<sub>i,j</sub>(1,1) 
</pre>

A function table consistent with the example preamble above would be:
<pre>
2
 -4.0 4.0
2
 -4.0 4.0
2 
 -4.0 4.0
2
 -4.0 4.0
4
# This is an additional comment
 1.5 -1.5 
 -1.5 1.5
4
 1.5 -1.5
 -1.5 1.5
4
 1.5 -1.5
 -1.5 1.5
</pre>
A full example file is given in test_ising.txt.
 