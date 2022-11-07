function x = fft2c(x,varargin)
%FFTC Centered discrete Fourier transform.
%	FFT2C(X) is the centered, discrete Fourier transform (DFT) of 
%	vector X.  For matrices, the FFT operation is applied to each
%	column. For N-D arrays, the FFT operation operates on the first
%	non-singleton dimension.
%	
%	FFT2C(X,N) is the N-point FFT, padded with zeros if X has less
%	than N points and truncated if it has more.
%	
%	FFT2C(X,[],DIM) or FFT2C(X,N,DIM) applies the FFT operation across
%	the dimension DIM.
%
%	See also IFFTC, FFT, IFFT, FFT2, IFFT2, FFTSHIFT.

%	v1.0, 6/25/01, Michael Butterfield, <gte881s@prism.gatech.edu>
%	v1.1, 6/26/01, Erik Zeek, <zeekec@mad.scientist.com>
%		Updated help.  Changed to VARARGIN implementation.
%
%	$Revision: 1.3 $ $Date: 2002/04/18 22:05:29 $
%
%	$Log: fftc.m,v $
%	Revision 1.3  2002/04/18 22:05:29  zeekec
%	Reused the input as the output.
%	
%	Revision 1.2  2001/07/10 01:10:00  zeekec
%	Library cleanup.  Added, deleted, and moved files.
%	
%

x = fftshift(fft2(ifftshift(x),varargin{:}));

