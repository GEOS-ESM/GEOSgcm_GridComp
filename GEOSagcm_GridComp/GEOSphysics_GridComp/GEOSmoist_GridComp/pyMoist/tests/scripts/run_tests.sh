# #!/bin/bash
# rm -rf ./.gt_cache_*
# export PACE_FLOAT_PRECISION=32
# export FV3_DACEMODE=Python
# python -m pytest -v -s -x\
#     --data_path=../../test_data/11.5.2/Moist/TBC_C24_L72_Debug \
#     --backend=dace:cpu \
#     --grid=default \
#     --multimodal_metric \
#     ..


#!/bin/bash
# rm -rf ./.gt_cache_*
export PACE_FLOAT_PRECISION=32
export FV3_DACEMODE=Python
python -m pytest -v -x \
    --data_path=/home/charleskrop/netcdfs \
    --backend=dace:cpu \
    --grid=default \
    --multimodal_metric \
    ..