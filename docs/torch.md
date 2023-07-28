>[!note]
>Always ensure to use `torch::kFloat64` as the `dtype` of any tensors, since libBigWig returns double precision (`double`) arrays. Use the provided `constants::tensor_opts` torch `TensorOpts` object in `proc_bigWigs.h` whenever possible during tensor initializations.

# Saving the Binned Tensors
[Loading a torch::tensor saved with C++, in Python](https://github.com/pytorch/pytorch/issues/39623#issuecomment-640662556)