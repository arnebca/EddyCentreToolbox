function data=nc_getall(filename)
% function to get all variables in a structure array using the new matlab NETCDF commands
% input string with datafile location and name
if nargout > 1
    disp(' nc_getall >> please use ONE or NO output file')
    return
end
try
    ncid = netcdf.open(filename,'nc_nowrite');
catch
    disp(' nc_getall >> unable to open input file')
    data = [];
    return
end

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for varid = 0:nvars-1
    try
        
        [varname,xtype,dimids,natts] = netcdf.inqVar(ncid,varid);
        
        if nargout == 0
            eval(['global ' varname]);
            if length(dimids) < 3
                if natts > 0
                    for attstep = 1:natts
                        attname = netcdf.inqAttName(ncid,varid,attstep-1);
                        if strcmp(attname,'scale_factor')
                            scale_factor = double(netcdf.getAtt(ncid,varid,'scale_factor'));
                        else
                            scale_factor = 1;
                        end
                        if strcmp(attname,'add_offset')
                            add_offset = double(netcdf.getAtt(ncid,varid,'add_offset'));
                        else
                            add_offset = 0;
                        end
                    end
                else
                    scale_factor = 1;
                    add_offset = 0;
                end
                evalin('base',[varname ' = double(netcdf.getVar(' num2str(ncid) ',' num2str(varid) ')'') .* ' num2str(scale_factor) ' + ' num2str(add_offset) ';']);
            else
                if natts > 0
                    for attstep = 1:natts
                        attname = netcdf.inqAttName(ncid,varid,attstep-1);
                        if strcmp(attname,'scale_factor')
                            scale_factor = double(netcdf.getAtt(ncid,varid,'scale_factor'));
                        else
                            scale_factor = 1;
                        end
                        if strcmp(attname,'add_offset')
                            add_offset = double(netcdf.getAtt(ncid,varid,'add_offset'));
                        else
                            add_offset = 0;
                        end
                    end
                else
                    scale_factor = 1;
                    add_offset = 0;
                end
                evalin('base',[varname ' = double(netcdf.getVar(' num2str(ncid) ',' num2str(varid) ')).* ' num2str(scale_factor) ' + ' num2str(add_offset) ';']);
            end
            
        else
            if length(dimids) < 3
                if natts > 0
                    for attstep = 1:natts
                        attname = netcdf.inqAttName(ncid,varid,attstep-1);
                        if strcmp(attname,'scale_factor')
                            scale_factor = double(netcdf.getAtt(ncid,varid,'scale_factor'));
                        else
                            scale_factor = 1;
                        end
                        if strcmp(attname,'add_offset')
                            add_offset = double(netcdf.getAtt(ncid,varid,'add_offset'));
                        else
                            add_offset = 0;
                        end
                    end
                else
                    scale_factor = 1;
                    add_offset = 0;
                end
                data.(varname) = double(netcdf.getVar(ncid,varid)'') .* scale_factor + add_offset;
                
            else
                if natts > 0
                    for attstep = 1:natts
                        attname = netcdf.inqAttName(ncid,varid,attstep-1);
                        if strcmp(attname,'scale_factor')
                            scale_factor = double(netcdf.getAtt(ncid,varid,'scale_factor'));
                        else
                            scale_factor = 1;
                        end
                        if strcmp(attname,'add_offset')
                            add_offset = double(netcdf.getAtt(ncid,varid,'add_offset'));
                        else
                            add_offset = 0;
                        end
                    end
                else
                    scale_factor = 1;
                    add_offset = 0;
                end
                data.(varname) = double(netcdf.getVar(ncid,varid)) .* scale_factor + add_offset;
                
            end
        end
    catch
        %disp([' nc_getall >> unable to get var: ' varname ])
    end
end
netcdf.close(ncid)



