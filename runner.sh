(echo "bcKink0 ic1 phi4kink LeapFrog"
 echo "0 0 0"               # bc parameters
 echo "0 5 2"               # ic parameters
 echo "5 1 0 0"             # v  parameters
 echo "2 -60 60"            # tf x0 xf
 echo "10000 5000")  |python3 main.py


sed -i '/^bcp1/,$!d' inc.py
cat inc.py