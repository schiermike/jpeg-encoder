Jpeg Encoder
============

This is a minimalistic implementation of a JPEG encoder for demonstration purposes. It has the following features/limitations:

 *  baseline-DCT profile
 *  fixed 4:2:0 chroma subsampling
 *  dynamic luma/chroma quantization
 *  only dimensions which are multiples of 16 are supported
 *  dynamic huffman table creation

 Only ``PPM`` image files are accepted as input.