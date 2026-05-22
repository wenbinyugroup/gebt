
# Inputs

GEBT has a build-in 1D finite element mesh generator, hence it is very easy to prepare the input file for a beam assembly.
For users to understand the meaning of the input data, we explain each data entry line by line here.

The first line lists three analysis control parameters arranged as:

```
analysis_flag  niter  nstep
```

where `analysis_flag` can be 0, 1, 2, 3.
If `analysis_flag` is 0, GEBT will carry out a static analysis.
If `analysis_flag` is 1, GEBT will carry out a steady state analysis (i.e., neglecting the time derivatives in the formulation).
If `analysis_flag` is 2, GEBT will carry out a time marching to obtain the transient time history of the dynamic response.
If `analysis_flag` is 3, GEBT will carry out an eigenvalue analysis to obtain the frequencies and mode shapes.
`niter` is the maximum number of iterations for the nonlinear analysis.
If `niter` is equal to 1, GEBT will compute results corresponding to a linear theory.
`nstep` is the number of steps.
For static analysis, it can be used to increase the load gradually to help the convergence.
For dynamic analysis, it is the number of time steps.

If `analysis_flag` is not equal to 0, the next two lines will be used to provide the angular velocity of the global frame a and the linear velocity of the starting pointing of the first member, and their corresponding time functions.
They are arranged as:

```
   w1     w2     w3
tf_w1  tf_w2  tf_w3
   v1     v2     v3
tf_v1  tf_v2  tf_v3
```

where `w1`, `w2`, `w3`, `v1`, `v2`, `v3` are real numbers denoting the magnitude, while `tf_w1`, `tf_w2`, `tf_w3`, `tf_v1`, `tf_v2`, `tf_v3` are integer numbers denoting the number of the corresponding time functions.
The values of these velocities are equal to the magnitude multiplied by values evaluated by the corresponding time function.
The unit of angular velocity will be rad/second.
If inch if chosen as the unit of length, the unit of linear velocity will be inch/second.

If `analysis_flag` is equal to 3, the next line will be used to provide an integer `nev` for the total number of frequencies and corresponding mode shapes to be extracted from the eigenvalue analysis.

The next line lists nine integers arranged as:

```
nkp  nmemb  ncond_pt  nmate  nframe  ncond_mb  ndistr  ntimefun  ncurv
```

where
- `nkp` is the total number of key points,
- `nmemb` the total number of members,
- `ncond_pt` the total number of points having prescribed conditions (including boundary conditions),
- `nmate` the total number of cross-sections,
- `nframe` the total number of frames,
- `ncond_mb` the total number of members having prescribed, distributed loadings,
- `ndistr` the total number of functions used to approximate the distributed loads,
- `ntimefun` the total number of time functions, `ncurv` total number of initial curvature/twist sets including $\boldsymbol { k } _ { 1 } , \boldsymbol { k } _ { 2 } , \boldsymbol { k } _ { 3 }$

The next nkp lines are the coordinates for each key point arranged as:

```
kp_no  x1  x2  x3
```

where `kp_no` is an integer representing the unique number assigned to each key point and `x1`, `x2`, `x3` are three real numbers describing the location $( x _ { 1 } , x _ { 2 } , x _ { 3 } )$ of the point.
Although the arrangement of `kp_no` is not necessary to be consecutive, every point starting from 1 to `nkp` should be present.

The next `nmemb` lines list eight integers for each member.
They are arranged as:

```
memb_no  kp_1  kp_2  mate_no1  mate_no2  frame_no  ndiv  curv_no
```

where `memb_no` is the number of member and `kp_1` is the starting point and `kp_2` is the ending point of the member. `mate_no1` is the cross-section number of the starting point of the member, `mate_no2` is the cross-section number of the ending point of the member.
If the member has a uniform cross-section, these two numbers will be the same.
Two different cross-sections can be assigned to one member, which implicitly assumes that the sectional properties are linearly varying along the length.
`frame_no` is the frame number of the starting point of the member.
If it is set to be 0, then the frame is the same as the global frame, i.e., $b_i = a_i$.
`ndiv` is the total number of elements you want to divide this member (the element is uniformly divided into `ndiv` segments), `curv_no` is the initial curvature/twist set number (0 means a straight member).
Although the arrangement of `memb_no` is not necessary to be consecutive, every member starting from 1 to `nmemb` should be present.

The next `ncond_pt` blocks define the point conditions of prescribed displacements, rotations, forces, or moments.
They are arranged as:

```
kp_no
dof_1  dof_2  dof_3  dof_4  dof_5  dof_6
val_1  val_2  val_3  val_4  val_5  val_6
 tf_1   tf_2   tf_3   tf_4   tf_5   tf_6
 ff_1   ff_2   ff_3   ff_4   ff_5   ff_6
```

where `kp_no` is the key point where the prescribed condition is applied, `dof_i` (`i` = 1, 2, ..., 6) are the prescribed degrees of freedom and they are six integers with values ranging from 1 to 12.
`val_i` (`i` = 1, 2, ..., 6) are six real numbers for the prescribed values for the corresponding degree of freedom.
`tf_i` (`i` = 1, 2, ..., 6) are six integer numbers for the corresponding number of the time function.
`ff_i` (`i` = 1, 2, ..., 6) are six integer numbers (could be either 0 or 1) indicating whether the corresponding prescribed quantity is a follower quantity or not.
If it is a follower, then the flag is set to be 1.
Otherwise, it is 0.
We assume that only forces/moments can be follower quantities and if there is a follower component at the point, GEBT assumes that the first three prescribed degrees of freedom can only be either 7, 8, 9 or 10, 11, 12.
If the assumption is too restrictive, please let the author know.
The prescribed value is calculated as val i multiplying the corresponding time function value.
For example if at point 2, we applied follower forces and dead moments with forces are constant with respect to time and moments are prescribed by time function 1 (defined later), then we have the following

```
 2
 7   8   9  10  11  12
F1  F2  F3  M1  M2  M3
 0   0   0   1   1   1
 1   1   1   0   0   0
```

where `F1`, `F2`, `F3`, `M1`, `M2`, `M3` are corresponding prescribed values, which could be zero, implying no external forces applied at this point.
In a beam assembly, the key points can be further classified as boundary points or connection points.
If a point is connected to only one member, it is a boundary point and six boundary conditions must be applied to this point.
If a point is connected to more than one member, it is a connection point.
The displacements/rotations of this point can be prescribed.
Concentrated forces/moments can also be applied to the connection point.

The next `nmate` blocks provide the structural and inertial properties for each cross-section.
First for the flexibility matrix, described using 36 real numbers arranged as:

```
mate_no
S_11  S_12  S_13  S_14  S_15  S_16
S_12  S_22  S_23  S_24  S_25  S_26
S_13  S_23  S_33  S_34  S_35  S_36
S_14  S_24  S_34  S_44  S_45  S_46
S_15  S_25  S_35  S_45  S_55  S_56
S_16  S_26  S_36  S_46  S_56  S_66
```

These values are defined in Eq. (3) and can be directly copied from VABS output file.
If one wants to input stiffness matrix instead, the user needs to use `preprocess.stiff.f90` to replace `preprocess.f90` and recompile the code.
If `analysis_flag` is not equal to 0, we also need to provide inertial properties represented by the mass matrix which is arranged as:

```
m_11  m_12  m_13  m_14  m_15  m_16
m_12  m_22  m_23  m_24  m_25  m_26
m_13  m_23  m_33  m_34  m_35  m_36
m_14  m_24  m_34  m_44  m_45  m_46
m_15  m_25  m_35  m_45  m_55  m_56
m_16  m_26  m_36  m_46  m_56  m_66
```

These values are defined in Eq. (4) and can be directly copied from VABS output file.

If `nframe` > 0, the next nframe blocks provide the direction cosine matrix for each member, which are described using nine real numbers arranged as:
```
frame_no
C11  C12  C13
C21  C22  C23
C31  C32  C33
```

These values are defined in Eq. (1).

If `ncond_mb` > 0, the next `ncond_mb` blocks define members having prescribed loads.
They are arranged as:

```
memb_no
 dn_1   dn_2   dn_3   dn_4   dn_5   dn_6
val_1  val_2  val_3  val_4  val_5  val_6
 tf_1   tf_2   tf_3   tf_4   tf_5   tf_6
 ff_1   ff_2   ff_3   ff_4   ff_5   ff_6
```

where
- `memb_no` is the member where the prescribed condition is applied,
- `dn_i` (`i` = 1, 2, ..., 6) are the corresponding number of distributed functions.
- `val_i` (`i` = 1, 2, ..., 6) are six real numbers for the prescribed values arranged corresponding to $f _ { 1 }$, $f _ { 2 }$, $f _ { 3 }$, $m _ { 1 }$, $m _ { 2 }$, $m _ { 3 }$ with $f _ { 1 }$, $f _ { 2 }$, $f _ { 3 }$ as three components of the distribution force and $m _ { 1 }$, $m _ { 2 }$, $m _ { 3 }$ as three moments of the distribution moment.
- `tf_i` (`i` = 1, 2, ..., 6) are six integer numbers for the corresponding number of the time function.
- `ff_i` (`i` = 1, 2, ..., 6) are six integer numbers (could be either 0 or 1) indicating whether the corresponding prescribed quantity is a follower quantity or not.
If it is a follower, then the flag is set to be 1.
Otherwise, it is 0.
The prescribed value is calculated as `val_i` multiplying the corresponding time function value and multiplying the corresponding distribution function value.

If `ndistr`> 0, the next `ndistr` blocks provide the distributed load functions.
They are arranged as:

```
fun_no
c0  c1  c2  c3  c4  c5
```

where `c_i` (`i` = 0, ..., 5) corresponds to the first five coefficients of the Chebychev polynomials used to approximate the given load function.
The function is calculated as $f(s) = \sum_{i=0}^{5} c_{i} T_{i}(s)$, where $s$ is the length along the member and $T _ { 0 } = 1$, $T _ { 1 } = s$, $T _ { 2 } = 2 s ^ { 2 } - 1$, $T _ { 3 } = 4 s ^ { 3 } - 3 s$, $T _ { 4 } = 8 s ^ { 4 } - 8 s ^ { 2 } + 1$, $T _ { 5 } = 16 s ^ { 5 } - 20 s ^ { 3 } + 5 s$.
For example, for a linearly distributed load $f ( s ) = 1000 + 500s$, we need to provide the input as `1000 500 0 0 0 0`.
Here GEBT assumes all given functions can be approximated by the Chebychev polynomials up to the fifth order.

If `ncurv`> 0, the next `ncurv` blocks provide the initial curvatures and twist.
They are arranged as:

```
curv_no
k1  k2  k3
```

where `k1` is the initial twist, `k2` is the initial curvature around $x _ { 2 }$ direction, and `k3` is the initial curvature around $x _ { 3 }$ direction.

If `ntimefun`> 0 or `analysis_flag`=2, the next line lists two numbers for the simulation range arranged as:

```
starting_time  ending_time
```

and the next `ntimefun` blocks define the time functions which are arranged as follows:

```
fun_no
fun_type
ts  te
[fun_block]
```

where
- `fun_no` is the time function number,
- `fun_type` is the time function type, and
- `ts`, `te` are the starting and ending time of function definition specifically.

Currently, user can define two types of time functions.
If `fun_type` is 0, it is a user defined as piecewise linear function based on functional values provided for different time instances.
The `[fun_block]` is arranged as follows:

```
n   
t1  f1   
t2  f2   
...   
tn  fn
```

where `n` is the number of entries needed for this time function.
The time entries $t _ { i }$ is arranged in an increasing fashion.
That is, we have $t _ { 2 } > t _ { 1 }$, $t _ { 3 } > t _ { 2 }$, and etc.
The time function is assumed to be piecewise linear defined by $t _ { i }$ and $f _ { i }$.
If $t < t _ { 1 }$, the function will remain the same as $f _ { 1 }$.
If $t > t _ { n }$, the function will remain the same as $f _ { n }$.
Note this type of time functions, simulation time, and load steps can be used to control the load increment in static analysis.

If `fun_type` is 1, the time function is defined as a summation of a series of harmonics as follows

$$
f ( t ) = \sum _ { i = 1 } ^ { N } h _ { i } ( t ) = \sum _ { i = 1 } ^ { N } a _ { i } \sin 2 \pi ( t / T _ { i } + \phi _ { i } )\tag{5}
$$

The `fun_block` is arranged as follows:

```
n
a_1  T_1  phi_1
a_2  T_2  phi_2
...
a_n  T_n  phi_n
```

where
- `a_i` is the magnitude,
- `T_i` is the period, and
- `phi_i` is the phase.

The input file should be ended with a blank line to avoid any possible incompatibility of different computer systems.
The input file can be given any name as long as the total number of the characters of the name including extension is not more than 256.
You are suggested to use a unique extension say get for you to identify such files with GEBT.
For the convenience of the user to identify mistakes in the input file, all the inputs are echoed in a file named `input_file_name.ech`.   
Error messages are also written at the end of `input_file_name.ech`.

If `analysis_flag`=2, a file named `input_file_name.ini` should also exist in the same directory.
If there are a total of `nelem` elements in the structure, this file contains 2 × `nelem` lines with the first `nelem` lines providing the initial positions and rotations $( u _ { 1 } \quad u _ { 2 } \quad u _ { 3 } \quad \theta _ { 1 } \quad \theta _ { 2 } \quad \theta _ { 3 } )$ for each element and the next `nelem` lines providing the corresponding derivatives $( { \dot { u } } _ { 1 } { \dot { u } } _ { 2 } { \dot { u } } _ { 3 } { \dot { \theta } } _ { 1 } { \dot { \theta } } _ { 2 } { \dot { \theta } } _ { 3 } )$ for each element.

