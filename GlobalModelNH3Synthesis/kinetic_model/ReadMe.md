# 0D kinetic model

This project generates a global kinetic model on the basis of an input file. 
This file must attain to the following rules, to ensure the working of the code;

Define a reaction as 

   `number,reaction formula,ReactionFunction(corresponding_constants)`
    
where, the number is merely for readability. 
The rules regarding the `reaction formula` are discussed in section 1.
The rules of the `ReactionFunction` are discussed in section 2.


## 1. Reaction Formula
For instance, 

      A + B -> C + D

* use `->` to differentiate between the reactants and products
* use `+` to add multiple reactants and/or products
* insert spaces around the aforementioned characters
* non-reacting species are not written down, i.e. electrons `e` and third party species `M`
* adsorbed species are denoted with `s`, e.g. `Ns`.


## 2. Reaction Function
...



## TODO
- include LookUp tables for the electron impact reactions
- ...


