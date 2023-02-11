#!/bin/bash

case "$1" in
    build)
        echo "Building 3PG-SoNWaL"
        package=`ls *.tar*`
        ## Command to install 3PG-SoNWaL package
        R --vanilla CMD build .
        ;;   
    install)
        echo "Installing 3PG-SoNWaL"
        package=`ls *.tar*`
        ## Command to install 3PG-SoNWaL package
        R --vanilla CMD INSTALL --no-lock $package
        ;;
    uninstall)
        echo "Removing 3PG-SoNWaL"
        ## Command for removing 3PG-SoNWaL package
        R --vanilla CMD REMOVE SoNWaL
        ;;
    *)
        echo "Usage ./sonwall.sh {build|install|uninstall}"
        exit 1
        ;;
esac

exit 0
        
