function [mov1,mov2,xShift,yShift] = correct_jitter_v2_2color(mov1,mov2,upScale,trunc)
% FUNCTION [MOV,XSHIFT,YSHIFT] = CORRECT_JITTER(MOV,UPSCALE,TRUNC)
%
% This is a function to quickly align movies with subpixel resolution
% to the first frame using a fft based approcah.  Starting with the second
% frame, it will attempt to calculate the xshift and yshift to align with
% the first frame, and then warp the second frame using those values.  It
% will then continue to align each following frame to the previously
% corrected version.
%
% Note this will only correct image translations, not shear or
% magnification issues.
%
% Created by: Joel Kralj - 2017/12/06
if nargin < 2
    upScale = 10;
    trunc = 2;
elseif nargin < 3
    trunc = 2;
end

numRegister = 5;

[ysize,xsize,nframes] = size(mov1);
xShift = zeros(nframes,1);
yShift = zeros(nframes,1);

yTrunc = round(ysize/2) - ...
    round(ysize/(2*trunc))+1:round(ysize/2)+round(ysize/(2*trunc))-1;
xTrunc = round(xsize/2) - ...
    round(xsize/(2*trunc))+1:round(xsize/2)+round(xsize/(2*trunc))-1;

% Start the alignment
alignImage = mov1(yTrunc,xTrunc,1);
alignRegister = zeros([size(alignImage) numRegister]);
alignRegister(:,:,1) = mov1(yTrunc,xTrunc,1);
for f = 2:nframes
    mov1(:,:,f) = imtranslate(mov1(:,:,f),[1.03*xShift(f-1) 1.03*yShift(f-1)]);
    mov2(:,:,f) = imtranslate(mov2(:,:,f),[1.03*xShift(f-1) 1.03*yShift(f-1)]);
    [yTmp,xTmp] = fftregister(fft2(alignImage),fft2(mov1(yTrunc,xTrunc,f)),upScale);
    mov1(:,:,f) = imtranslate(mov1(:,:,f),[1.0*xTmp 1.0*yTmp]);
    mov2(:,:,f) = imtranslate(mov2(:,:,f),[1.0*xTmp 1.0*yTmp]);
    xShift(f) = xTmp+xShift(f-1);
    yShift(f) = yTmp+yShift(f-1);
    for g = numRegister-1:1
        alignRegister(:,:,g+1) = alignRegister(:,:,g);
    end
    alignRegister(:,:,1) = mov1(yTrunc,xTrunc,f);
    alignImage = mean(alignRegister,3);
end
end

function [row_shift,col_shift] = fftregister(buf1ft,buf2ft,usfac)
if ~exist('usfac','var')
    usfac = 1;
end

[nr,nc]=size(buf2ft);
Nr = ifftshift(-fix(nr/2):ceil(nr/2)-1);
Nc = ifftshift(-fix(nc/2):ceil(nc/2)-1);


% Start with usfac == 2
CC = ifft2(FTpad(buf1ft.*conj(buf2ft),[2*nr,2*nc]));
CCabs = abs(CC);
[row_shift, col_shift] = find(CCabs == max(CCabs(:)),1,'first');
CCmax = CC(row_shift,col_shift)*nr*nc;
% Now change shifts so that they represent relative shifts and not indices
Nr2 = ifftshift(-fix(nr):ceil(nr)-1);
Nc2 = ifftshift(-fix(nc):ceil(nc)-1);
row_shift = Nr2(row_shift)/2;
col_shift = Nc2(col_shift)/2;
% If upsampling > 2, then refine estimate with matrix multiply DFT
if usfac > 2
    % Initial shift estimate in upsampled grid
    row_shift = round(row_shift*usfac)/usfac;
    col_shift = round(col_shift*usfac)/usfac;
    dftshift = fix(ceil(usfac*1.5)/2); %% Center of output array at dftshift+1
    % Matrix multiply DFT around the current shift estimate
    CC = conj(dftups(buf2ft.*conj(buf1ft),ceil(usfac*1.5),ceil(usfac*1.5),usfac,...
        dftshift-row_shift*usfac,dftshift-col_shift*usfac));
    % Locate maximum and map back to original pixel grid
    CCabs = abs(CC);
    [rloc, cloc] = find(CCabs == max(CCabs(:)),1,'first');
    CCmax = CC(rloc,cloc);
    rloc = rloc - dftshift - 1;
    cloc = cloc - dftshift - 1;
    row_shift = row_shift + rloc/usfac;
    col_shift = col_shift + cloc/usfac;
end

if nr == 1
    row_shift = 0;
end
if nc == 1
    col_shift = 0;
end

end

function out=dftups(in,nor,noc,usfac,roff,coff)
% function out=dftups(in,nor,noc,usfac,roff,coff);
% Upsampled DFT by matrix multiplies, can compute an upsampled DFT in just
% a small region.
% usfac         Upsampling factor (default usfac = 1)
% [nor,noc]     Number of pixels in the output upsampled DFT, in
%               units of upsampled pixels (default = size(in))
% roff, coff    Row and column offsets, allow to shift the output array to
%               a region of interest on the DFT (default = 0)
% Recieves DC in upper left corner, image center must be in (1,1) 
% Manuel Guizar - Dec 13, 2007
% Modified from dftus, by J.R. Fienup 7/31/06

% This code is intended to provide the same result as if the following
% operations were performed
%   - Embed the array "in" in an array that is usfac times larger in each
%     dimension. ifftshift to bring the center of the image to (1,1).
%   - Take the FFT of the larger array
%   - Extract an [nor, noc] region of the result. Starting with the 
%     [roff+1 coff+1] element.

[nr,nc]=size(in);
% Set defaults
if exist('roff', 'var')~=1, roff=0;  end
if exist('coff', 'var')~=1, coff=0;  end
if exist('usfac','var')~=1, usfac=1; end
if exist('noc',  'var')~=1, noc=nc;  end
if exist('nor',  'var')~=1, nor=nr;  end
% Compute kernels and obtain DFT by matrix products
kernc=exp((-1i*2*pi/(nc*usfac))*( ifftshift(0:nc-1).' - floor(nc/2) )*( (0:noc-1) - coff ));
kernr=exp((-1i*2*pi/(nr*usfac))*( (0:nor-1).' - roff )*( ifftshift([0:nr-1]) - floor(nr/2)  ));
out=kernr*in*kernc;
end

function [ imFTout ] = FTpad(imFT,outsize)
% imFTout = FTpad(imFT,outsize)
% Pads or crops the Fourier transform to the desired ouput size. Taking 
% care that the zero frequency is put in the correct place for the output
% for subsequent FT or IFT. Can be used for Fourier transform based
% interpolation, i.e. dirichlet kernel interpolation. 
%
%   Inputs
% imFT      - Input complex array with DC in [1,1]
% outsize   - Output size of array [ny nx] 
%
%   Outputs
% imout   - Output complex image with DC in [1,1]
% Manuel Guizar - 2014.06.02

if ~ismatrix(imFT)
    error('Maximum number of array dimensions is 2')
end
Nout = outsize;
Nin = size(imFT);
imFT = fftshift(imFT);
center = floor(size(imFT)/2)+1;

imFTout = zeros(outsize);
centerout = floor(size(imFTout)/2)+1;

% imout(centerout(1)+[1:Nin(1)]-center(1),centerout(2)+[1:Nin(2)]-center(2)) ...
%     = imFT;
cenout_cen = centerout - center;
imFTout(max(cenout_cen(1)+1,1):min(cenout_cen(1)+Nin(1),Nout(1)),max(cenout_cen(2)+1,1):min(cenout_cen(2)+Nin(2),Nout(2))) ...
    = imFT(max(-cenout_cen(1)+1,1):min(-cenout_cen(1)+Nout(1),Nin(1)),max(-cenout_cen(2)+1,1):min(-cenout_cen(2)+Nout(2),Nin(2)));

imFTout = ifftshift(imFTout)*Nout(1)*Nout(2)/(Nin(1)*Nin(2));
end
