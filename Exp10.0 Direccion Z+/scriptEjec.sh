#!/bin/sh

if [[ $1 == 'listo' ]];then
    echo 'segunda ejecucion aqui va a terminar xq no se hace mas'
    #time ./exampleB1 Tarea.mac # ya estoy usando Main.mac
    time ./exampleB1 Main.mac
    #./exampleB1
else
    echo 'primera ejecucion'
    cmake ..
    make -j3
    chmod +x  scriptEjec.sh
    bash scriptEjec.sh 'listo'
fi


