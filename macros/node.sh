#!bin/bash
root -l debug7_test_grain.root <<EOC
tree1->MakeClass("myNode")
.q
EOC
rm myNode.C
