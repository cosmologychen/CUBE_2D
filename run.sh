# 重新编译、生成 IC、再运行主程序与后处理
make \
  && ./utilities/ic.x > ic.log \
  && ./main.x > main.log \
  && ./utilities/density_power.x  > density_power.log \
  && ./utilities/dsp.x  > dsp.log \
  && ./utilities/lightcone.x > lightcone.log


tail -n 30 main.log
ls -lh   /mnt/18T/output_2D/lc_c10/3000_1024_onerun/0.000_delta_c.bin