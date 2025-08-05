# orthotropic_pf_subroutine
This subroutine implements a phase-field (PF) model with multiple damage variables in Abaqus using a combination of UMAT and UEL subroutines.
Below are the steps to correctly set up the Abaqus model.

To configure the model, you can create a standard Abaqus .inp file following conventional procedures. Only a few specific details need to be taken into account.

First, the material must be defined as a User Material, and four parameters must be specified: E1, E2, nu12 and G12. In addition, the Depvar option must be included, specifying 12 state variables.
```inp
*Material, name=MechMaterial
*User Material, constants=4
E1, E2, nu12, G12
*Depvar
12
```

Second, the step must be of type Coupled Temperature–Displacement, and the Steady State option should be selected.
```inp
*Step, name=Step-1, inc=nInc
*Coupled Temperature-displacement, creep=none, steady state
initialTimeInc, totalTimeInc, minimumTimeInc, maximumTimeInc
```

The mesh must be generated in Abaqus, taking into account that the subroutine is designed to work only with linear quadrilateral elements. Once this is done and the boundary conditions have been applied, the .inp file can be created.

In the .inp file, the number of elements must be doubled. This is done by expanding the element list, adding as many new elements as there are in the original model. It is important to ensure that if the original model contains n elements, then element i + n must have the same connectivity as element i, for all i.

The first n elements must be assigned the user element type. This element type must be defined in the .inp file and then assigned to these elements accordingly. The assigned degrees of freedom in the user element definition correspond to the temperature degrees of freedom in Abaqus.
```inp
*User element, nodes=4, type=U1, properties=9, coordinates=2, variables=1
 11,12,13,14,15,16
*Element, type=U1
```

Starting from element n + 1, an Abaqus element must be assigned. Specifically, the CPS4 element.
```inp
*Element, type=CPS4
```

Create an element set for the phase-field (PF) element layer (to solve the damage problem) and another element set for the Abaqus element layer (to solve the mechanical problem).
```inp
*Elset, elset=SetPF, generate
 1,  n,  1
*Elset, elset=SetMech, generate
 n+1,  2*n,  1
```

The properties to be used with the user element in the subroutine must be assigned. There are 9 properties to define: E1, E2, nu12, G12, l01, Gc1, l02, Gc2, and theta. The angle theta corresponds to the angle formed by direction 1 with respect to the global X-axis and should be given in degrees.
```inp
*UEL PROPERTY, ELSET=SetPF
E1, E2, nu12, G12, l01, Gc1, l02, Gc2
theta
```

Also, make sure to assign the User Material properties to the element set used for the mechanical problem, along with the material orientation theta.
```inp
*Orientation, name=Ori-1
1., 0., 0., 0., 1., 0.
3, theta
** Section: Section-1
*Solid Section, elset=SetMech, orientation=Ori-1, material=MechMaterial
,
```

Finally, it is necessary to modify the value of the variable nelem in the module kvisual of the subroutine so that it matches the number of elements in one material layer, i.e., n.
