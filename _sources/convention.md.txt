
# Conventions

To understand the inputs and interpret outputs of the program correctly, we need to explain some conventions used by GEBT.

Firstly, GEBT uses a right-hand coordinate system, the beam coordinate system denoted as $x _ { 1 } , x _ { 2 }$ and $x _ { 3 } .$ , where $x _ { 1 }$ is along the beam axis and $x _ { 2 }$ and $x _ { 3 }$ are the local Cartesian coordinates of the cross section, see Figure 1 for a beam with an arbitrary cross section.
Unit vectors $\mathbf { b } _ { i }$ are the corresponding base vectors.
The cross-sectional coordinates $( x _ { 2 }$ and $x _ { 3 } )$ are only needed for the purpose to define its orientation.
It is noted that the beam coordinate system remains the same as the undeformed beam coordinate system defined in Ref. [5] for straight beam members.
For curved member, the undeformed beam coordinate system is changing along the beam axis and the user needs to provide the undeformed beam coordinate system at the starting point of the beam member and the changing undeformed beam coordinate system along with the beam axis will be calculated automatically by the code.
Usually, for a beam assembly, we define a global coordinate system, say with the triad ${ \bf a } _ { i } ,$ then the orientation of a beam member can relate to the global coordinate system using a direction cosine matrix

$$
{ \left\{ \begin{array} { l } { \mathbf { a } _ { 1 } } \\ { \mathbf { a } _ { 2 } } \\ { \mathbf { a } _ { 3 } } \end{array} \right\} } = { \left[ \begin{array} { l l l } { C _ { 1 1 } } & { C _ { 1 2 } } & { C _ { 1 3 } } \\ { C _ { 2 1 } } & { C _ { 2 2 } } & { C _ { 2 3 } } \\ { C _ { 3 1 } } & { C _ { 3 2 } } & { C _ { 3 3 } } \end{array} \right] } \left\{ \begin{array} { l } { \mathbf { b } _ { 1 } } \\ { \mathbf { b } _ { 2 } } \\ { \mathbf { b } _ { 3 } } \end{array} \right\} = C ^ { a b } \left\{ \begin{array} { l } { \mathbf { b } _ { 1 } } \\ { \mathbf { b } _ { 2 } } \\ { \mathbf { b } _ { 3 } } \end{array} \right\}
$$

```{figure} images/vabs_beam_csys.png

VABS beam coordinate system.
```

Secondly, GEBT uses two points to specify each member with $x _ { 1 }$ of the member pointing to the second point.
For dynamic analysis, there are 18 values for each element within a member arranged as

$$
\left[ { \begin{array} { l l l l l l l l l l l l l } { u _ { 1 } } & { u _ { 2 } } & { u _ { 3 } } & { \theta _ { 1 } } & { \theta _ { 2 } } & { \theta _ { 3 } } & { F _ { 1 } } & { F _ { 2 } } & { F _ { 3 } } & { M _ { 1 } } & { M _ { 2 } } & { M _ { 3 } } & { P _ { 1 } } & { P _ { 2 } } & { P _ { 3 } } & { H _ { 1 } } & { H _ { 2 } } & { H _ { 3 } } \end{array} } \right] ^ { T }
$$

Thirdly, GEBT assumes that the cross-sectional properties have already been calculated by VABS and given in the form of a flexibility matrices and mass matrix.
The flexibility matrix is defined according to the following equation

$$
\left\{ \begin{array} { c } { { \gamma _ { 1 1 } } } \\ { { 2 \gamma _ { 1 2 } } } \\ { { 2 \gamma _ { 1 3 } } } \\ { { \kappa _ { 1 } } } \\ { { \kappa _ { 2 } } } \\ { { \kappa _ { 3 } } } \end{array} \right\} = \left[ \begin{array} { c c c c c c } { { S _ { 1 1 } } } & { { S _ { 1 2 } } } & { { S _ { 1 3 } } } & { { S _ { 1 4 } } } & { { S _ { 1 5 } } } & { { S _ { 1 6 } } } \\ { { S _ { 1 2 } } } & { { S _ { 2 2 } } } & { { S _ { 2 3 } } } & { { S _ { 2 4 } } } & { { S _ { 2 5 } } } & { { S _ { 2 6 } } } \\ { { S _ { 1 3 } } } & { { S _ { 2 3 } } } & { { S _ { 3 3 } } } & { { S _ { 3 4 } } } & { { S _ { 3 5 } } } & { { S _ { 3 6 } } } \\ { { S _ { 1 4 } } } & { { S _ { 2 4 } } } & { { S _ { 3 4 } } } & { { S _ { 4 4 } } } & { { S _ { 4 5 } } } & { { S _ { 4 6 } } } \\ { { S _ { 1 5 } } } & { { S _ { 2 5 } } } & { { S _ { 3 5 } } } & { { S _ { 4 5 } } } & { { S _ { 5 5 } } } & { { S _ { 5 6 } } } \\ { { S _ { 1 6 } } } & { { S _ { 2 6 } } } & { { S _ { 3 6 } } } & { { S _ { 4 6 } } } & { { S _ { 5 6 } } } & { { S _ { 6 6 } } } \end{array} \right] \left[ \begin{array} { c } { { F _ { 1 } } } \\ { { F _ { 2 } } } \\ { { F _ { 3 } } } \\ { { M _ { 1 } } } \\ { { M _ { 2 } } } \\ { { M _ { 3 } } } \end{array} \right]
$$

where $\gamma _ { 1 1 }$ is the beam axial stretching strain measure, $2 \gamma _ { 1 2 }$ and $2 \gamma _ { 1 3 }$ are the engineering transverse shear strains along $x _ { 2 }$ and $x _ { 3 }$ respectively, $\kappa _ { 1 }$ is the twist measure, and $\kappa _ { 2 }$ and $\kappa _ { 3 }$ are the curvature measures around $x _ { 2 }$ and $x _ { 3 }$ respectively.
Note that any of the values in the flexibility matrix $( S _ { i j } )$ can be zero to eliminate any of the deformation mechanism which the user believes to be negligible.
For example, to analyze a pure bending around $x _ { 2 }$ , one can set $S _ { 5 5 }$ to be the bending flexibility and other terms to be zeroes.

The elements of the mass matrix are arranged as

$$
\left[ \begin{array} { c c c c c c } { \mu } & { 0 } & { 0 } & { 0 } & { \mu x _ { m 3 } } & { - \mu x _ { m 2 } } \\ { 0 } & { \mu } & { 0 } & { - \mu x _ { m 3 } } & { 0 } & { 0 } \\ { 0 } & { 0 } & { \mu } & { \mu x _ { m 2 } } & { 0 } & { 0 } \\ { 0 } & { - \mu x _ { m 3 } } & { \mu x _ { m 2 } } & { i _ { 2 2 } + i _ { 3 3 } } & { 0 } & { 0 } \\ { \mu x _ { m 3 } } & { 0 } & { 0 } & { 0 } & { i _ { 2 2 } } & { - i _ { 2 3 } } \\ { - \mu x _ { m 2 } } & { 0 } & { 0 } & { 0 } & { - i _ { 2 3 } } & { i _ { 3 3 } } \end{array} \right]
$$

where $\mu$ is mass per unit length, $( x _ { m 2 } , x _ { m 3 } )$ is the location of mass center, $i _ { 2 2 }$ is the mass moment of inertia about $x _ { 2 }$ axis, i33 is the mass moment of inertia about $x _ { 3 }$ axis, i23 is the product of inertia.

GEBT allows users to use various kinds of units.
However, it is necessary to be absolutely consistent in the choice of units to avoid errors.
Particularly, users must never use the pound as a unit of mass to avoid confusion.
When pounds are used for force and feet for length, the unit of mass must be the s ${ \mathrm { l u g } } = { \mathrm { l b } } { \mathrm { - s e c } } ^ { 2 } / { \mathrm { f t } } ;$ if inches are used for length along with pounds for force, then the unit of mass must be $\scriptstyle \mathrm { l b - s e c ^ { 2 } / i n }$

