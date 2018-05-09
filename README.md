# FEM_code
Some code about FEM.

I have improved some algorithms in the process of learning. So the code in different examples may use different algorithms to solve the same problems. And usually the relative complex element examples would show some efficient algorithms.

1 Improvement in assembling global stiffness matrix. Use the index to assemble efficiently.

2 Improvement in solving unknown variables by rearranging and spliting matrix.

3 Using the application of parametric elements to simplify some complex interpolation and integration.

-------------Updated 3/20/2018-----------------------

Most of the codes were refined recently. I have changed the input and output mode. First using interactive mode, but then I was inspired that reading input file and yielding output file could be better.

There were also some trivial changes in using 3-dimensional matrices, solving global stiffness equations by partioning matrix, etc.

The most important one is to formulate the procedure of making isoparametric element models, both in 2d and 3d form. Actually some of the steps in formulating finite elements are still not so clear in my mind. I have to add comments in Chinese.

I'm sure that all of the codes have to be improved incessantly.

--------------Updated 5/9/2018------------------------

I have put all in codes wriiten by MATLAB into another folder in this repository. And I began to use Python to rewrite most of the codes. Still, there will be some improvements in this edition.
