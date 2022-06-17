g++ -fopenmp \
-I /opt/X11/include -L /opt/X11/lib -lX11 \
main.cpp -o main \
&& rm -f \
data/con/*.vtk \
data/phi/*.vtk \
data/con/*.csv \
data/phi/*.csv \
figures/con_xy/*.png \
figures/con_yz/*.png \
figures/con/*.png \
figures/phi/*.png \
&& ./main \
&& rm main \
&& python plot1d.py