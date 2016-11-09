;	PLOT_NTW_GR_RELATION.PRO	
;
;	Plots Gutenberg-Richter relation for the catalog  
;
;	
;	July 30, 2010
;
;	+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;
    pro plot_ntw_gr_relation,                                                       $
        mag_ge, number_of_mag_bins, maximum_magnitude, year_file, mag_file,         $
		Time_End, Time_Start, plot_time_step, total_time_steps,                     $
		mag_max, mag_min, mag_step, number_of_small_eq_array, fname, kount
;
;	***************************************************************
;
;
;	Now plot the results --  
;
;	Load color table
	loadct, 39	; rainbow + white
;
	!P.POSITION =  [.15, .15, .95, .90]
	
;
;	FIRST PLOT: Try pdf first
;


  	!P.charsize = 1.0

    calc_gr_arrays,                                             $
        mag_file, mag_min, mag_step,number_of_mag_bins,         $
        total_time_steps, number_of_small_eq_array,             $
        pdf_array, hst_array, cdf_gr_array, gr_array, mag_array

    plot_array   =   dblarr(number_of_mag_bins)
    plot_array(*) = 0.0

    plot_array(*) = gr_array(*)

;
;   *******************************************************
;
;     Fit the line to find the b-value
;

    xr    =   fltarr(number_of_mag_bins)
    yr    =   fltarr(number_of_mag_bins)

      print,  ' '
      print,  ' Number of points is',number_of_mag_bins
      print,  ' '
      print,  ' Enter LOW and HIGH mag. values'
      print,  ' to determine b-value (magnitudes)'
      print,  ' (to bypass this put HIGH < LOW)'

      mag1 = 0.
      mag2 = 0.
      read,  mag1, mag2

      k = -1
       for i = 0,number_of_mag_bins-1 do begin
       if(mag_array(i) ge mag1 and mag_array(i) le mag2) then begin
       k = k + 1
       xr(k) = mag_array(i)
       yr(k) = plot_array(i)
       print, ' k, xr, yr: ', k, xr(k), yr(k)
       endif

       endfor
      kmax = k+1

      if(kmax gt 1) then begin

      linfit, xr, yr, slope, cept, errs, errc, kmax,s2

      bval = - slope
      print,  ' '
      print,  ' b - value is:'
      print,   bval,'   +/-',errs
      print, ' '
      print, ' '

      endif

;
;   ********************************************************
;
    n_ge_4 = 0
    n_ge_5 = 0
    n_ge_6 = 0
    n_ge_7 = 0
    n_ge_8 = 0

        for i=0L, kount-1L do begin
        if (mag_file(i) ge 4.0) then n_ge_4=n_ge_4 + 1
        if (mag_file(i) ge 5.0) then n_ge_5=n_ge_5 + 1
        if (mag_file(i) ge 6.0) then n_ge_6=n_ge_6 + 1
        if (mag_file(i) ge 7.0) then n_ge_7=n_ge_7 + 1
        endfor

        print, ' '
        print, 'Total Number EQs > 4: ', n_ge_4
        print, 'Total Number EQs > 5: ', n_ge_5
        print, 'Total Number EQs > 6: ', n_ge_6
        print, 'Total Number EQs > 7: ', n_ge_7
        print, ' '
    
;
;   ********************************************************
;

    print, 'mag_array: ', mag_array(*)
    print, 'pdf_array: ', pdf_array(*)
    print, 'gr_array: ', gr_array(*)
    print, 'hst_array: ', hst_array(*)

    Max_Plot = float(long(max(plot_array))) + 1.0

    Max_Plot = 3.0
	Min_Plot = 0.0

    X_min = 5.0
    X_max = 8.0

;   X_min = 2.0
;   X_max = 6.0
    
;    X_max = float(long(max(mag_file))) + 1.0


;    X_max = 7.0

     plot,[0],[0], /ynoz, $                          ; dummies
        xstyle=1,ystyle=1, $              ; force exact axis range
        xrange=[X_min,X_max], $   
        yrange=[Min_Plot, Max_Plot], $
		title   = 'Gutenberg-Richter Relation', $
		xtitle  = 'Magnitude', $
		ytitle  = 'Log10 (GR Frequency)', charsize = 1.25


;
;	***********************************************************************************************
;
;	Make filled circles
;
        A = findgen(16) * (!pi*2/16.)
        usersym, cos(A), sin(A), /fill

		oplot,mag_array,plot_array,psym=8,symsize=0.85,color=50


;   Plot lines
;
        solid_line_x=fltarr(2)
        solid_line_y=fltarr(2)

        dashed_line_x=fltarr(2)
        dashed_line_y=fltarr(2)

        slope=(yr(kmax-1)-yr(0))/(xr(kmax-1)-xr(0))
        b_intercept = yr(kmax-1) - slope*xr(kmax-1)

        solid_line_x(0) = xr(0)
        solid_line_y(0) = yr(0)

        solid_line_x(1) = xr(kmax-1)
        solid_line_y(1) = yr(kmax-1)

        dashed_line_x(0) = xr(kmax-1)
        dashed_line_y(0) = yr(kmax-1)

        dashed_line_x(0) = x_min
        dashed_line_y(0) = slope*x_min + b_intercept

        dashed_line_x(1) = - b_intercept/slope
        dashed_line_y(1) = 0.0

            if (mag1 lt mag2) then begin
        	    oplot,solid_line_x,solid_line_y,linestyle=0,color=250,thick = 3.0
        	    oplot,dashed_line_x,dashed_line_y,linestyle=2,color=250,thick = 3.0
            	endif
        
;

;   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;
;
	legends = 'on'

	if(legends eq 'on') then begin

        xs = 0.05*(X_max - X_min)  +  X_min
        ys = 0.05*(Max_Plot - Min_Plot) + Min_Plot
        astr3 = 'Catalog: '
        astr11 = strcompress(fname)
        astr12 = astr3 + astr11
        xyouts, xs, ys, astr12, charsize = 0.5

        xs = 0.05*(X_max - X_min)  +  X_min
        ys = 0.115*(Max_Plot - Min_Plot) + Min_Plot
        astr1 = 'Years: '
        astr2 = strcompress(string(Time_Start))
        astr21 = strmid(astr2,1,7)
		astr3 = strcompress(string(Time_End))
		astr31 = strmid(astr3,1,7)
		astr4 = astr1 + astr21 + '-' + astr31		
        xyouts, xs, ys, astr4, charsize = 0.85

      if(kmax gt 1) then begin

        xs = 0.05*(X_max - X_min)  +  X_min
        ys = 0.18*(Max_Plot - Min_Plot) + Min_Plot
        astr1 = 'b-Value: '
        astr2 = strcompress(string(bval))
        astr21 = strmid(astr2,1,5)
		astr3 = strcompress(string(errs))
		astr31 = strmid(astr3,1,5)
		astr4 = astr1 + astr21 + '  +/- ' + astr31		
        xyouts, xs, ys, astr4, charsize = 1.0

      endif


	endif
;
;	&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&



;
; opportunity to stop while in this environment
;
 
   dummy = '' & read,'Press Return to continue...  ',dummy
;
;
;	End of code
;
	end


	
	
