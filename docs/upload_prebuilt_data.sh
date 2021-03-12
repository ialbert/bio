set -uex

REMOTE=www@biostarhandbook.com:/home/www/book/data_www/bio

rsync $FLAGS ~/.bio/taxdb.json  $HOST:$DIR
rsync $FLAGS ~/.bio/taxdb.sqlite  $HOST:$DIR