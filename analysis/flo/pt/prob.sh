:> test_im1h.dat; for i in {1..59}; do grep "Prob.*IM1H" ../out/nvt_$i.out | cut -d":" -f2 >> test_im1h.dat; done
:> test_hoac.dat; for i in {1..59}; do grep "Prob.*HOAC" ../out/nvt_$i.out | cut -d":" -f2 >> test_hoac.dat; done
