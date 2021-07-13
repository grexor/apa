if [ -d `pwd`/data ]
then
  chmod 777 data
else
  mkdir data
  chmod 777 data
fi

docker run -v `pwd`/data:/home/apauser/data --user apauser --hostname groovy -ti apa bash --login
