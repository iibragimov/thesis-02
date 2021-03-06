!  Console1.f90 
!****************************************************************************
!
!  PROGRAM: ��������� ��������� � ������ ���������. 
!
!  PURPOSE: ������� ��� ������ ����� ��������� ������.
!
!  DESCRIPTION: ����� ����� ��-��� ��������, �������� ����� ��������� �� ��������� �� �������������.
!
!  NOTICE: 
!           1. ����� v(s) �������� �� ������� ��������, ��������� ����� ��
!
!           2. ������, ��� �������� � �� ��������� ���������, ������ ��� �� ��������� ��������� gamma o
!
!           3. ��������� �� ����� ����������, � ����������� ���������� 4, � �� ����� �������� � ����� �������
!
!           4. 
!
!           5. ��������� ������� ds_dgamma
!  
!  DONE:
!           1. ������� ��� �����_�, ��� � � |C|      
!
!           2. ������� ��� �_�, �������� �� �����  
!
!           3. ������� ��� ������� � ����������
!
!
!****************************************************************************

    program Console1     
    !����������� ������ � ����������� � �������������� �������������� (IMSL)
    use mod
    
    !������ ��� ������ � ����
    port = 1
    
    ! ��������� �����������
    ga = 7.0d0/6.0d0*pi
    gb = pi/3.0d0
    gc = 3.3d0/d2*pi
    gm = 7.5d0/4.0d0*pi
    mod_C = d5
    
    !����� ��������    
    call find_const(ga, gb, gc, gm, mod_C)
    
    !������������� �������� �� ��������
    call v_s
    
    !����� ����
    call current_lines
    
    ! ������ � ���� ���������� ��������
    call find_shape_of_plast
    
    end program
    
   

