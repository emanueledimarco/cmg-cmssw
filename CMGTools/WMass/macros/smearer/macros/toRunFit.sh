#!/bin/zsh    

writeScript() {   

     root -b<<EOF

.L fitEleResponseMap.C++ 
fitEleResponseMap()  
.q
EOF     
EOF

}

writeScript;