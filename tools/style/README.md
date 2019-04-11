# Coding convention

The style for source code follows the [LLVM
standard](https://llvm.org/docs/CodingStandards.html).


## This coding style

There are few deviations from the above coding style, defined in the file
`format` in this directory.

### Variable naming convention

The [naming conventions](https://llvm.org/docs/CodingStandards.html#name-types-functions-variables-and-enumerators-properly)
follow the LLVM standard with one exception:

* **Variable names** They are not treated as 'Nouns'.  Consequently, we do not
write them with capitalized starting letters, nor do we write them in
CamelCase.  Long variable names may use underscores for better readability.
Examples:
```
int mycount;
int mycount_reverse;
float pressure_drag; // - for a long scope where the whole code may not fit on screen.
float pdrag;         // - for a short scope that fits on +/-entire screen.  Such a
                     // variable declaration must be accompanied by a comment.
float pd;            // - same rule as for 'pdrag' above.  According to the LLVM
                     // standard, 'pdrag' should be preferred.
```
