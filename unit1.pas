unit Unit1;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, FileUtil, Forms, Controls, Graphics, Dialogs, ExtCtrls,
  StdCtrls, Grids, Buttons, ComCtrls;

type

  { TForm1 }

  TForm1 = class(TForm)
    Bevel1: TBevel;
    BTN_waves_clear: TButton;
    BTN_signal_xy_clear: TButton;
    BTN_signal_xy_template1: TButton;
    BTN_signal_xy_template2: TButton;
    BTN_spectrum_syntesis: TButton;
    BTN_waves_template1: TButton;
    BTN_spectrum_analisys: TButton;
    BTN_waves_N: TButton;
    Btn_reset_t: TButton;
    BTN_waves_rnd: TButton;
    BTN_syntes_img_waves: TButton;
    BTN_signal_Nx_Ny: TButton;
    BTN_waves_template2: TButton;
    CB_timer: TCheckBox;
    CB_waves_legend: TCheckBox;
    CB_z_level: TCheckBox;
    CB_ab_phi_model: TCheckBox;
    CB_z_waves: TCheckBox;
    Edit_signal_Ny: TEdit;
    Edit_waves_width: TEdit;
    Edit_waves_sx: TEdit;
    Edit_waves_height: TEdit;
    Edit_waves_sy: TEdit;
    Edit_w1: TEdit;
    Edit_t: TEdit;
    Edit_dt: TEdit;
    Edit_Waves_N: TEdit;
    Edit_k_dist: TEdit;
    Edit_signal_Nx: TEdit;
    Label1: TLabel;
    Label10: TLabel;
    Label12: TLabel;
    LBL_img_waves_x: TLabel;
    Label11: TLabel;
    LBL_img_waves_y: TLabel;
    LBL_img_waves_XY_val: TLabel;
    Label2: TLabel;
    Label3: TLabel;
    Label4: TLabel;
    Label5: TLabel;
    Label6: TLabel;
    Label7: TLabel;
    Label8: TLabel;
    Label9: TLabel;
    PBox: TPaintBox;
    SG_waves: TStringGrid;
    SG_signal_xy: TStringGrid;
    Timer1: TTimer;
    UpDown_sy: TUpDown;
    UpDown_sx: TUpDown;
    procedure BTN_waves_clearClick(Sender: TObject);
    procedure Btn_reset_tClick(Sender: TObject);
    procedure BTN_syntes_img_wavesClick(Sender: TObject);
    procedure BTN_waves_rndClick(Sender: TObject);
    procedure BTN_signal_Nx_NyClick(Sender: TObject);
    procedure BTN_signal_xy_template1Click(Sender: TObject);
    procedure BTN_signal_xy_template2Click(Sender: TObject);
    procedure BTN_spectrum_analisysClick(Sender: TObject);
    procedure BTN_signal_xy_clearClick(Sender: TObject);
    procedure BTN_spectrum_syntesisClick(Sender: TObject);
    procedure BTN_waves_NClick(Sender: TObject);
    procedure BTN_waves_template1Click(Sender: TObject);
    procedure BTN_waves_template2Click(Sender: TObject);
    procedure CB_ab_phi_modelChange(Sender: TObject);
    procedure CB_timerChange(Sender: TObject);
    procedure CB_waves_legendChange(Sender: TObject);
    procedure CB_z_levelChange(Sender: TObject);
    procedure CB_z_wavesChange(Sender: TObject);
    procedure FormCreate(Sender: TObject);
    procedure PBoxMouseMove(Sender: TObject; Shift: TShiftState; X, Y: Integer);
    procedure PBoxPaint(Sender: TObject);
    procedure Timer1Timer(Sender: TObject);
    procedure UpDown_sxClick(Sender: TObject; Button: TUDBtnType);
    procedure UpDown_syClick(Sender: TObject; Button: TUDBtnType);
  private
    procedure Draw_img_waves;
    procedure print_params;
    procedure get_params;
  public

  end;

type
TComplex=record
   re,im:real;
end;
type TComplexVector=array of TComplex;
type TComplexMatrix=array of TComplexVector;

type
  t_wave=record
    xb,yb:real;
    x0,y0:real;
    fx,fy:real;
    amp:real;
    phi:real;
    a,b:real;
  end;

var
  Form1: TForm1;
  init_flag:boolean;
  legend_flag, Z_level_flag, Z_waves_flag, ab_flag:boolean;
  MouseX, MouseY: integer; MouseShift:TShiftState;
  SummaryBitmap:TBitmap;
  img_width,img_height:integer;

  t,dt:real;
  waves_sx,waves_sy:real;
  waves_width,waves_height:integer;
  k_dist,dx,dy:real;
  N_waves: integer;
  waves: array[0..10000] of t_wave;
  IMG_waves: array[0..1023,0..1023]of real;

  signal_nx, signal_ny: integer;
  signal_xy:array[0..1023,0..1023]of real;
  fourie_w1:real;

implementation

{$R *.lfm}

function c_add(a,b:TComplex):TComplex;
begin c_add.re:=a.re+b.re; c_add.im:=a.im+b.im; end;

function c_mul(a,b:TComplex):TComplex;
begin c_mul.re:=a.re*b.re-a.im*b.im; c_mul.im:=a.re*b.im+a.im*b.re; end;

function c_root_of_one_CW(k,n:integer):TComplex;
var phi:real;
begin
     phi:=-2*PI*k/n;
     c_root_of_one_CW.re:=cos(phi);
     c_root_of_one_CW.im:=sin(phi);
end;

function c_root_of_one_CCW(k,n:integer):TComplex;
var phi:real;
begin
     phi:=2*PI*k/n;
     c_root_of_one_CCW.re:=cos(phi);
     c_root_of_one_CCW.im:=sin(phi);
end;

//преобразование Фурье (спектральный анализ)
//медленный алгоритм со степенью сложности O(N*N)
//память для DFT_t и DFT_f должна быть выделена до вызова подпрограммы
procedure DFT_analysis_slow(var DFT_t,DFT_f:TComplexVector);
var k,j,N,M:integer;
begin
  N:=length(DFT_t);
  M:=length(DFT_f);
  for k:=0 to M-1 do
  begin
     DFT_f[k].re:=0; DFT_f[k].im:=0;
     for j:=0 to N-1 do
        DFT_f[k]:=c_add(DFT_f[k],c_mul(DFT_t[j],c_root_of_one_CW(k*j,N)));
     DFT_f[k].re:=DFT_f[k].re/M;
     DFT_f[k].im:=DFT_f[k].im/M;
  end;
end;

//преобразование Фурье (спектральный синтез)
//медленный алгоритм со степенью сложности O(N*N)
//память для DFT_t и DFT_f должна быть выделена до вызова подпрограммы
procedure DFT_syntesis_slow(var DFT_f,DFT_t:TComplexVector);
var k,j,N,M:integer;
begin
  N:=length(DFT_t);
  M:=length(DFT_f);
  for j:=0 to N-1 do
  begin
    DFT_t[j].re:=0; DFT_t[j].im:=0;
    for k:=0 to M-1 do
      DFT_t[j]:=c_add(DFT_t[j],c_mul(DFT_f[k],c_root_of_one_CCW(k*j,N)));
  end;
end;

//двумерное преобразование Фурье (спектральный анализ) - медленный способ
procedure DFT_analysis_2D_slow(var DFT_t,DFT_f:TComplexMatrix);
var k1,k2,j1,j2,img_height,img_width:integer;
begin
  img_height:=length(DFT_t);
  if img_height>0 then
  begin
    img_width:=length(DFT_t[0]);
    if img_width>0 then
    for k1:=0 to img_height-1 do
    for k2:=0 to img_width-1 do
    begin
      DFT_f[k1,k2].re:=0; DFT_f[k1,k2].im:=0;
      for j1:=0 to img_height-1 do
      for j2:=0 to img_width-1 do
        DFT_f[k1,k2]:=c_add(DFT_f[k1,k2],
                            c_mul(DFT_t[j1,j2],
                              c_mul(c_root_of_one_CW(k1*j1,img_height),
                                c_root_of_one_CW(k2*j2,img_width)))
                            );
      DFT_f[k1,k2].re:=DFT_f[k1,k2].re/(img_width*img_height);
      DFT_f[k1,k2].im:=DFT_f[k1,k2].im/(img_width*img_height);
    end;
  end;
end;

//двумерное преобразование Фурье (спектральный синтез) - медленный способ
procedure DFT_syntesis_2D_slow(var DFT_f,DFT_t:TComplexMatrix);
var k1,k2,j1,j2,img_height,img_width:integer;
begin
  img_height:=length(DFT_f);
  if img_height>0 then
  begin
    img_width:=length(DFT_f[0]);
    if img_width>0 then
    for j1:=0 to img_height-1 do
    for j2:=0 to img_width-1 do
    begin
      DFT_t[j1,j2].re:=0; DFT_t[j1,j2].im:=0;
      for k1:=0 to img_height-1 do
      for k2:=0 to img_width-1 do
        DFT_t[j1,j2]:=c_add(DFT_t[j1,j2],
                            c_mul(DFT_f[k1,k2],
                              c_mul(c_root_of_one_CCW(k1*j1,img_height),
                                c_root_of_one_CCW(k2*j2,img_width)))
                            );
    end;
  end;
end;
//========================================================================

procedure calc_ab_from_phiamp(var w:t_wave);
begin
    w.a:=w.amp*cos(w.phi);
    w.b:=-w.amp*sin(w.phi);
end;

procedure calc_phiamp_from_ab(var w:t_wave);
begin
    w.amp:=sqrt(sqr(w.a)+sqr(w.b));
    if w.a<>0
    then
      begin
        w.phi:=abs(arctan(w.b/w.a));
        if (w.a<0)and(w.b<0) then w.phi:=pi-w.phi;
        if (w.a>0)and(w.b>0) then w.phi:=-w.phi;
        if (w.a<0)and(w.b>0) then w.phi:=-pi+w.phi;
      end
    else
      begin
        if w.b>0 then w.phi:=pi/2;
        if w.b<0 then w.phi:=-pi/2;
        if w.b=0 then w.phi:=0;
      end;
end;

procedure correction_ab_phi;
var k:integer;
begin
  if ab_flag then for k:=0 to N_waves-1 do calc_phiamp_from_ab(waves[k])
             else for k:=0 to N_waves-1 do calc_ab_from_phiamp(waves[k]);
end;

procedure clr_waves;
var k:integer;
begin
  for k:=0 to N_waves-1 do
  begin
    waves[k].xb:=0;
    waves[k].yb:=0;
    waves[k].x0:=0;
    waves[k].y0:=0;
    waves[k].amp:=0;
    waves[k].fx:=0;
    waves[k].fy:=0;
    waves[k].phi:=0;
    waves[k].a:=0;
    waves[k].b:=0;
  end;
end;

procedure clr_signal;
var i,j:integer;
begin
  for i:=0 to signal_nx-1 do
  for j:=0 to signal_ny-1 do
    signal_xy[i,j]:=0;
end;

procedure spectrum_analisys;
var i,j,k:integer; DFT_f,DFT_t:TComplexMatrix;
begin
  setlength(DFT_f,signal_nx,signal_ny);
  setlength(DFT_t,signal_nx,signal_ny);

  for i:=0 to signal_nx-1 do
  for j:=0 to signal_ny-1 do
  begin
    DFT_t[i,j].re:=signal_xy[i,j];
    DFT_t[i,j].im:=0;
  end;

  DFT_analysis_2D_slow(DFT_t,DFT_f);

//  for k:=0 to N_waves-1 do
  for i:=0 to signal_nx-1 do
  for j:=0 to signal_ny-1 do
  begin
    k:=j*signal_nx+i;
    waves[k].fx:=fourie_w1*i;
    waves[k].fy:=fourie_w1*j;
    waves[k].a:=DFT_f[i,j].re;
    waves[k].b:=DFT_f[i,j].im;
    calc_phiamp_from_ab(waves[k]);
  end;

  setlength(DFT_f,0);
  setlength(DFT_t,0);
end;

procedure spectrum_syntesis;
var i,j,k:integer; DFT_f,DFT_t:TComplexMatrix; w:t_wave;
begin
  setlength(DFT_f,signal_nx,signal_ny);
  setlength(DFT_t,signal_nx,signal_ny);

//  for k:=0 to N_waves-1 do
  for i:=0 to signal_nx-1 do
  for j:=0 to signal_ny-1 do
  begin
    k:=j*signal_nx+i;
    w:=waves[k];
    if not(ab_flag) then calc_ab_from_phiamp(w);
    DFT_f[i,j].re:=w.a;
    DFT_f[i,j].im:=w.b;
  end;

  DFT_syntesis_2D_slow(DFT_f,DFT_t);

  for i:=0 to signal_nx-1 do
  for j:=0 to signal_ny-1 do
    signal_xy[i,j]:=DFT_t[i,j].re;

  setlength(DFT_f,0);
  setlength(DFT_t,0);
end;

function wave_elliptic(w:t_wave; tx,ty:real):real;
begin wave_elliptic:=w.a*cos(tx*2*pi*w.fx+w.phi)*cos(ty*2*pi*w.fy+w.phi);   end;

procedure gen_img_waves;
var k,x,y:integer; w:t_wave; tx,ty,tmp_amp,xr,yr:real;
begin
  for x:=0 to waves_width-1 do
  for y:=0 to waves_height-1 do
  begin
    IMG_waves[x,y]:=0;
    xr:=x/k_dist+waves_sx/k_dist;
    yr:=y/k_dist+waves_sy/k_dist;
    for k:=0 to N_waves-1 do
    begin
      tx:=xr-waves[k].xb;
      ty:=yr-waves[k].yb;
      tmp_amp:=wave_elliptic(waves[k],t+tx,t+ty);
      IMG_waves[x,y]:=IMG_waves[x,y]+tmp_amp;
    end;
  end;
end;

{ TForm1 }
procedure TForm1.get_params;
var i,j,k:integer;
begin
  legend_flag:=CB_waves_legend.Checked;
  z_level_flag:=CB_z_level.Checked;
  ab_flag:=CB_ab_phi_model.Checked;
  z_waves_flag:=CB_z_waves.Checked;

  N_waves:=StrToInt(Edit_Waves_N.text);
  fourie_w1:=StrToFloat(Edit_w1.text);

  waves_sx:=StrToFloat(Edit_waves_sx.text);
  waves_sy:=StrToFloat(Edit_waves_sy.text);

  waves_width:=StrToInt(Edit_waves_width.text);
  waves_height:=StrToInt(Edit_waves_height.text);
  dx:=img_width div waves_width;
  dy:=img_height div waves_height;

  t:=StrToFloat(Edit_t.text);
  dt:=StrToFloat(Edit_dt.text);

  k_dist:=StrToFloat(Edit_k_dist.text);

  for k:=0 to N_waves-1 do
  begin
    waves[k].xb:=StrToFloat(SG_waves.Cells[1,k+1]);
    waves[k].yb:=StrToFloat(SG_waves.Cells[2,k+1]);
    waves[k].fx:=StrToFloat(SG_waves.Cells[3,k+1]);
    waves[k].fy:=StrToFloat(SG_waves.Cells[4,k+1]);
    waves[k].amp:=StrToFloat(SG_waves.Cells[5,k+1]);
    waves[k].phi:=StrToFloat(SG_waves.Cells[6,k+1]);
    waves[k].a:=StrToFloat(SG_waves.Cells[7,k+1]);
    waves[k].b:=StrToFloat(SG_waves.Cells[8,k+1]);
  end;
  correction_ab_phi;

  signal_nx:=StrToInt(Edit_signal_Nx.text);
  signal_ny:=StrToInt(Edit_signal_Ny.text);

  for i:=0 to signal_nx-1 do
  for j:=0 to signal_ny-1 do
    signal_xy[i,j]:=StrToFloat(SG_signal_xy.Cells[i+1,j+1]);
end;

procedure TForm1.print_params;
var i,j,k:integer;
begin
  CB_waves_legend.Checked:=legend_flag;
  CB_z_level.Checked:=z_level_flag;
  CB_ab_phi_model.Checked:=ab_flag;
  CB_z_waves.Checked:=z_waves_flag;

  Edit_Waves_N.text:=IntToStr(N_waves);

  Edit_waves_sx.text:=FloatToStr(waves_sx);
  Edit_waves_sy.text:=FloatToStr(waves_sy);

  Edit_waves_width.text:=IntToStr(waves_width);
  Edit_waves_height.text:=IntToStr(waves_height);

  Edit_t.text:=FloatToStrF(t,ffFixed,1,5);
  Edit_dt.text:=FloatToStr(dt);
  Edit_k_dist.text:=FloatToStrF(k_dist,ffFixed,1,3);
  Edit_w1.text:=FloatToStrF(fourie_w1,ffFixed,1,3);

  for k:=0 to N_waves-1 do
  begin
    SG_waves.Cells[0,k+1]:=IntToStr(k);
    SG_waves.Cells[1,k+1]:=FloatToStrF(waves[k].xb,ffFixed,1,5);
    SG_waves.Cells[2,k+1]:=FloatToStrF(waves[k].yb,ffFixed,1,5);
    SG_waves.Cells[3,k+1]:=FloatToStrF(waves[k].fx,ffFixed,1,5);
    SG_waves.Cells[4,k+1]:=FloatToStrF(waves[k].fy,ffFixed,1,5);
    SG_waves.Cells[5,k+1]:=FloatToStrF(waves[k].amp,ffFixed,1,5);
    SG_waves.Cells[6,k+1]:=FloatToStrF(waves[k].phi,ffFixed,1,5);
    SG_waves.Cells[7,k+1]:=FloatToStrF(waves[k].a,ffFixed,1,5);
    SG_waves.Cells[8,k+1]:=FloatToStrF(waves[k].b,ffFixed,1,5);
    SG_waves.Cells[9,k+1]:=FloatToStrF(waves[k].x0,ffFixed,1,5);
    SG_waves.Cells[10,k+1]:=FloatToStrF(waves[k].y0,ffFixed,1,5);
  end;

  Edit_signal_Nx.text:=IntToStr(signal_nx);
  Edit_signal_Ny.text:=IntToStr(signal_ny);

  for i:=0 to signal_nx-1 do SG_signal_xy.Cells[i+1,0]:=IntToStr(i);
  for i:=0 to signal_ny-1 do SG_signal_xy.Cells[0,i+1]:=IntToStr(i);

  for i:=0 to signal_nx-1 do
  for j:=0 to signal_ny-1 do
    SG_signal_xy.Cells[i+1,j+1]:=FloatToStrF(signal_xy[i,j],ffFixed,1,1);
end;

procedure TForm1.Draw_img_waves;
var k,x,y:integer; c:byte; tmp:LongInt; dst_ptr:PByte; dst_bpp:integer;
  amp,tx,ty,xr,yr:real;
begin
  SummaryBitmap.BeginUpdate(false);
  dst_bpp:=SummaryBitmap.RawImage.Description.BitsPerPixel div 8;
  dst_ptr:=SummaryBitmap.RawImage.Data;
  for y:=img_height-1 downto 0 do
  for x:=0 to img_width-1 do
  begin
    C:=0; tmp:=trunc(IMG_waves[trunc(x/dx),trunc(y/dy)])+128;
    if (tmp>=0)and(tmp<=255) then C:=tmp;
    if tmp<0 then C:=0;
    if tmp>255 then C:=255;
    dst_ptr^:=C; (dst_ptr+1)^:=C; (dst_ptr+2)^:=C;
    inc(dst_ptr,dst_bpp);
  end;
  SummaryBitmap.EndUpdate(false);
  PBox.Canvas.Draw(0,0,SummaryBitmap);

  if legend_flag then
  begin
    PBox.Canvas.Pen.Color:=clRed;
    PBox.Canvas.Brush.Style:=bsClear;
    PBox.Canvas.Font.Color:=clRed;
    for k:=0 to N_waves-1 do
    begin
      x:=trunc(dx*k_dist*(waves[k].xb-waves_sx/k_dist));
      y:=trunc(dy*k_dist*(waves[k].yb-waves_sy/k_dist));
      PBox.Canvas.EllipseC(x,PBox.Height-y,trunc(waves[k].amp),trunc(waves[k].amp));
      PBox.Canvas.TextOut(x,PBox.Height-y,IntToStr(k));
    end;
  end;

  if z_level_flag or z_waves_flag then
  begin
    PBox.Canvas.Pen.Color:=clYellow;
    for k:=-10 to 10 do
    begin
      PBox.Canvas.TextOut(0,PBox.Height div 2 - k*25,IntToStr(k*25));
      PBox.Canvas.Line(0,PBox.Height div 2 - k*25,PBox.Width-1,PBox.Height div 2 - k*25);
    end;
    PBox.Canvas.Pen.Color:=clRed;
    PBox.Canvas.Line(0,MouseY,PBox.Width-1,MouseY);
  end;

  if z_level_flag then
  begin
    PBox.Canvas.Pen.Color:=clBlue;
    PBox.Canvas.MoveTo(0,PBox.Height div 2);
    for x:=0 to waves_width-1 do
    begin
      amp:=IMG_waves[x,trunc((PBox.Height-MouseY)/dy)];
      PBox.Canvas.LineTo(trunc(dx*x),PBox.Height div 2 - trunc(amp));
    end;
  end;

  if z_waves_flag then
  begin
    yr:=(PBox.height-MouseY)/(dy*k_dist)+waves_sy/k_dist;
    for k:=0 to N_waves-1 do
    begin
      PBox.Canvas.Pen.Color:=clGreen;
      PBox.Canvas.MoveTo(0,PBox.Height div 2);
      ty:=yr-waves[k].yb;
      for x:=0 to waves_width-1 do
      begin
        xr:=x/k_dist+waves_sx/k_dist;
        tx:=xr-waves[k].xb;
        amp:=wave_elliptic(waves[k],t+tx,t+ty);
        PBox.Canvas.LineTo(trunc(x*dx),PBox.Height div 2 - trunc(amp));
      end;
    end;
  end;
end;

procedure TForm1.BTN_spectrum_analisysClick(Sender: TObject);
begin
  get_params;
  spectrum_analisys;
  print_params;
  get_params;
  print_params;
  gen_img_waves;
  Draw_img_waves;
end;

procedure TForm1.BTN_spectrum_syntesisClick(Sender: TObject);
begin
  get_params;
  spectrum_syntesis;
  print_params;
  get_params;
  print_params;
end;

procedure TForm1.BTN_syntes_img_wavesClick(Sender: TObject);
begin
  get_params;
  gen_img_waves;
  Draw_img_waves;
  print_params;
end;

procedure TForm1.Btn_reset_tClick(Sender: TObject);
begin
  t:=0; Edit_t.text:=FloatToStr(t);
  get_params;
  print_params;
  gen_img_waves;
  Draw_img_waves;
end;

procedure TForm1.BTN_waves_NClick(Sender: TObject);
begin
  N_waves:=StrToInt(Edit_Waves_N.text);
  SG_waves.RowCount:=N_waves+1;
  print_params;
  get_params;
  print_params;
  gen_img_waves;
  Draw_img_waves;
end;

procedure TForm1.BTN_waves_template1Click(Sender: TObject);
var k:integer;
begin
  clr_waves;
  for k:=0 to 0 do
  begin
    waves[k].xb:=0;
    waves[k].yb:=0;
    waves[k].fx:=1;
    waves[k].fy:=1;
    waves[k].amp:=100;
    waves[k].phi:=0;
    calc_ab_from_phiamp(waves[k]);
  end;
  print_params;
  get_params;
  print_params;
  gen_img_waves;
  Draw_img_waves;
end;

procedure TForm1.BTN_waves_template2Click(Sender: TObject);
var k:integer;
begin
  clr_waves;
  for k:=0 to N_waves-1 do
  begin
    waves[k].xb:=(random(waves_width)-(waves_width div 2))/k_dist;
    waves[k].yb:=0;
    waves[k].fx:=random(5)+1;
    waves[k].fy:=random(5)+1;
    waves[k].amp:=random(20)+20;
    waves[k].phi:=random()*2*pi-pi;
    calc_ab_from_phiamp(waves[k]);
  end;
  print_params;
  get_params;
  print_params;
  gen_img_waves;
  Draw_img_waves;
end;

procedure TForm1.BTN_waves_clearClick(Sender: TObject);
begin
  clr_waves;
  print_params;
  get_params;
  print_params;
  gen_img_waves;
  Draw_img_waves;
end;

procedure TForm1.BTN_waves_rndClick(Sender: TObject);
var k:integer;
begin
  for k:=0 to N_waves-1 do
  begin
    waves[k].xb:=(random()*waves_width-waves_width/2)/k_dist;
    waves[k].yb:=(random()*waves_height-waves_height/2)/k_dist;
    waves[k].fx:=random(5)+1;
    waves[k].fy:=random(5)+1;
    waves[k].amp:=50*random();
    waves[k].phi:=2*pi*random()-pi;
    calc_ab_from_phiamp(waves[k]);
  end;
  print_params;
  get_params;
  print_params;
  gen_img_waves;
  Draw_img_waves;
end;

procedure TForm1.BTN_signal_Nx_NyClick(Sender: TObject);
begin
  signal_nx:=StrToInt(Edit_signal_Nx.text);
  signal_ny:=StrToInt(Edit_signal_Ny.text);
  SG_signal_xy.ColCount:=signal_nx+1;
  SG_signal_xy.RowCount:=signal_ny+1;
  print_params;
end;

procedure TForm1.BTN_signal_xy_clearClick(Sender: TObject);
begin
  clr_signal;
  print_params;
end;

procedure TForm1.BTN_signal_xy_template1Click(Sender: TObject);
var i:integer;
begin
  clr_signal;
  for i:=0 to signal_nx-1 do signal_xy[i,0]:=50;
  print_params;
end;

procedure TForm1.BTN_signal_xy_template2Click(Sender: TObject);
var i,j:integer;
begin
  clr_signal;

  for i:=2 to signal_nx-3 do
  for j:=2 to signal_ny-3 do
    signal_xy[i,j]:=trunc(50*random()+50);

  for i:=3 to signal_nx-4 do
  for j:=3 to signal_ny-4 do
    signal_xy[i,j]:=0;

  print_params;
end;

procedure TForm1.CB_ab_phi_modelChange(Sender: TObject);
begin
  get_params;
  print_params;
  gen_img_waves;
  Draw_img_waves;
end;

procedure TForm1.CB_timerChange(Sender: TObject);
begin
  if CB_timer.Checked then Timer1.Enabled:=true else Timer1.Enabled:=false;
end;

procedure TForm1.CB_waves_legendChange(Sender: TObject);
begin
  get_params;
  gen_img_waves;
  Draw_img_waves;
end;

procedure TForm1.CB_z_levelChange(Sender: TObject);
begin
  MouseX:=255; MouseY:=255;
  get_params;
  gen_img_waves;
  Draw_img_waves;
end;

procedure TForm1.CB_z_wavesChange(Sender: TObject);
begin
  MouseX:=255; MouseY:=255;
  get_params;
  gen_img_waves;
  Draw_img_waves;
end;

procedure TForm1.FormCreate(Sender: TObject);
begin
  img_width:=PBox.width;
  img_height:=PBox.height;
  SummaryBitmap:=TBitmap.Create;
  SummaryBitmap.SetSize(img_width,img_height);

  timer1.Enabled:=false;
  N_waves:=1024;
  signal_nx:=32;
  signal_ny:=32;
  fourie_w1:=1;
  k_dist:=64;
  dt:=-0.1;

  waves_width:=128;
  waves_height:=128;

  waves_sx:=0;//-waves_width/2;
  waves_sy:=0;//-waves_height/2;

  dx:=img_width/waves_width;
  dy:=img_height/waves_height;

  SG_waves.RowCount:=N_waves+1;
  SG_waves.Cells[0,0]:='';
  SG_waves.Cells[1,0]:='Xb';
  SG_waves.Cells[2,0]:='Yb';
  SG_waves.Cells[3,0]:='fx';
  SG_waves.Cells[4,0]:='fy';
  SG_waves.Cells[5,0]:='amp';
  SG_waves.Cells[6,0]:='phi';
  SG_waves.Cells[7,0]:='a';
  SG_waves.Cells[8,0]:='b';
  SG_waves.Cells[9,0]:='x0';
  SG_waves.Cells[10,0]:='y0';
  clr_waves;

  SG_signal_xy.ColCount:=signal_nx+1;
  SG_signal_xy.RowCount:=signal_ny+1;

  SG_signal_xy.Cells[0,0]:='';
  clr_signal;

  print_params;
  get_params;
  print_params;
  gen_img_waves;
  //Draw_img_waves;
  init_flag:=true;
end;

procedure TForm1.PBoxMouseMove(Sender: TObject; Shift: TShiftState; X,Y: Integer);
var xr,yr:real;
begin
  MouseX:=X; MouseY:=Y; MouseShift:=Shift;
  xr:=MouseX/(dx*k_dist)+waves_sx/k_dist;
  yr:=(PBox.Height-MouseY)/(dy*k_dist)+waves_sy/k_dist;
  LBL_img_waves_x.caption:='x='+FloatToStrF(xr,ffFixed,1,1);
  LBL_img_waves_y.caption:='y='+FloatToStrF(yr,ffFixed,1,1);
  LBL_img_waves_XY_val.caption:='amp='+
                FloatToStrF(IMG_waves[trunc(MouseX/dx),
                                      trunc((PBox.Height-MouseY)/dy)],
                            ffFixed,1,2);

  if Z_level_flag or Z_waves_flag then Draw_img_waves;
end;

procedure TForm1.PBoxPaint(Sender: TObject);
begin
  if init_flag then
  begin
    gen_img_waves;
    Draw_img_waves;
  end;
end;

procedure TForm1.Timer1Timer(Sender: TObject);
begin
  t:=t+dt; Edit_t.Text:=floatToStr(t);
  get_params;
  gen_img_waves;
  if not(CB_z_level.checked) then Draw_img_waves;
  PBoxMouseMove(self,MouseShift,MouseX,MouseY);
  print_params;
end;

procedure TForm1.UpDown_sxClick(Sender: TObject; Button: TUDBtnType);
begin
  Edit_waves_sx.Text:=IntToStr(UpDown_sx.Position);
  get_params;
  gen_img_waves;
  Draw_img_waves;
  print_params;
end;

procedure TForm1.UpDown_syClick(Sender: TObject; Button: TUDBtnType);
begin
  Edit_waves_sy.Text:=IntToStr(UpDown_sy.Position);
  get_params;
  gen_img_waves;
  Draw_img_waves;
  print_params;
end;

end.

