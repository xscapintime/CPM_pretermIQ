rsync -avp -e 'ssh -p 28903' --include='*.csv' --include='*.npy' --exclude='*' 'liyang@172.18.4.20:/export/home/liyang/cpm/pt/*' .

