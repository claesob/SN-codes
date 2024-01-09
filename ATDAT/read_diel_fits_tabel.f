      character*120 ch
      open(11,file='diel_fits_badnell_orig.dat',status='old')
      open(6,file='diel_fits_badnell.dat')      
      do j=1,10000
         read(11,99)ch
c         write(6,99)ch
 99      format(a120)
         if(j>=3) then
c               write(6,*)j,i,ch(i:i)
            if((ch(1:3).eq.'  6'.or.ch(1:3).eq.'  7'.or.ch(1:3).eq.
     &           '  8'.or.ch(1:3).eq.' 10'.or.ch(1:3).eq.' 11'.or.
     &           ch(1:3).eq.' 12'.or.ch(1:3).eq.' 13'.or.ch(1:3).
     &           eq.' 14'.or.ch(1:3).eq.' 16'.or.ch(1:3).eq.' 18'.
     &           or.ch(1:3).eq.' 20'.or.ch(1:3).eq.' 26'.or.ch(1:3).
     &           eq.'  2').and.ch(9:9).eq.'1') then
               do i=120,1,-1                  
                  if(ch(i:i).ne.' ') then
                     num=i
                     if(num==34) then
                        n=7
                        write(6,91)n,ch(1:num), 0.0000, 0.0000, 0.0000 
     &                       ,0.0000,0.0000, 0.0000, 0.0000
 91                     format(i3,a34,1pe11.3,10e11.3)
                        goto 11
                     elseif(num==45) then
                        n=6
                        write(6,9)n,ch(1:num+1), 0.0000, 0.0000, 0.0000, 
     &                       0.0000,0.0000, 0.0000
 9                      format(i3,a46,1pe11.3,10e11.3)
                        goto 11
                     elseif(num==56) then
                        n=5
                        write(6,92)n,ch(1:num+1), 0.0000, 0.0000, 0.0000 
     &                       ,0.0000,0.0000
 92                     format(i3,a57,1pe11.3,10e11.3)
                        goto 11
                     elseif(num==67) then
                        n=4
                        write(6,93)n,ch(1:num+1), 0.0000, 0.0000, 0.0000 
     &                       ,0.0000
 93                     format(i3,a68,1pe11.3,10e11.3)
                        goto 11
                     elseif(num==78) then
                        n=3
                        write(6,94)n,ch(1:num+1), 0.0000, 0.0000, 0.0000
 94                     format(i3,a79,1pe11.3,10e11.3)
                        goto 11
                     elseif(num==89) then
                        n=2
                        write(6,95)n,ch(1:num+1), 0.0000, 0.0000 
 95                     format(i3,a90,1pe11.3,10e11.3)
                        goto 11
                     elseif(num==100) then
                        n=1
                        write(6,96)n,ch(1:num+1), 0.0000 
 96                     format(i3,a101,1pe11.3,10e11.3)
                        goto 11
                     elseif(num==111) then
                        n=0
                        write(6,97)n,ch(1:num+1)
 97                     format(i3,a112,1pe11.3,10e11.3)
                        goto 11                                                               
                     endif
                  endif
               enddo
            endif   
         endif         
 11      continue
      enddo
      end
      
