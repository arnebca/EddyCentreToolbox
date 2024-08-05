function varargout = lon360to180(l360,varargin)
% 
% [lon,mat1,mat2,...] = lon360to180(l360,'sort',mat1,mat2,...)

l180 = l360;
l180(l360 > 180) = l360(l360 > 180) - 360;

if nargin > 1
    
    if strcmpi(varargin{1} , 'sort') && nargin > 2
        
        [l180 , si] = sort(l180);
        varargout{1} = l180;
        
        % loop over additional input matrices
        for ni = 2 : nargin-1
            
            mat = varargin{ni};

            londim = find(size(mat) == length(l360));
            if length(londim) > 1
                warning('The additional input matrix has the same length as the longitude vector in more than one dimension. Make sure that the the matrix is re-ordered in the right way!');
                londim = find(size(l180) > 1);
                assert(size(l180,londim) == size(mat,londim) , 'Error: Not clear which is the longitude dimension of the additional input matrix! To solve this problem make sure that the longitude dimension of the matrix is the same as the dimension/orientation of the longitude vector.');
            end
            assert(~isempty(londim),'Error: The additional input matrix does not have the same length as the longitude vector in any dimension!');

            
            switch londim
                case 1
                    mat = mat(si,:,:,:);
                case 2
                    mat = mat(:,si,:,:);
                case 3
                    mat = mat(:,:,si,:);
                case 4
                    mat = mat(:,:,:,si);
            end
            
            varargout{ni} = mat;
            
        end
        
    elseif strcmp(varargin{1} , 'sort')
        
        l180 = sort(l180);
        varargout{1} = l180;
        
    else
        
    end
   
else
    
    varargout{1} = l180;
    
end
    
end