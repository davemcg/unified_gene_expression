- https://www.digitalocean.com/community/tutorials/how-to-set-up-shiny-server-on-ubuntu-14-04
- https://www.rstudio.com/products/shiny/download-server/

Commands run on Cyclops:
```
sudo su - -c "R -e \"install.packages('shiny', repos='http://cran.rstudio.com/')\""
sudo apt-get install gdebi-core
wget https://download3.rstudio.org/ubuntu-12.04/x86_64/shiny-server-1.4.2.786-amd64.deb
sudo gdebi shiny-server-1.4.2.786-amd64.deb
sudo ufw allow 3838 #open up shiny port
sudo su - -c "R -e \"install.packages('rmarkdown', repos='http://cran.rstudio.com/')\""
```


