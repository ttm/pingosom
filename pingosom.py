#-*- coding: utf-8 -*-
import numpy as n, random, os, sys, time
from scipy.io import wavfile as w
tfoo=time.time()
H=n.hstack
V=n.vstack

f_a = 44100. # Hz, frequência de amostragem

############## 2.2.1 Tabela de busca (LUT)
Lambda_tilde=Lt=1024.*16

# Senoide
foo=n.linspace(0,2*n.pi,Lt,endpoint=False)
S_i=n.sin(foo) # um período da senóide com T amostras

# Quadrada:
Q_i=n.hstack(  ( n.ones(Lt/2)*-1 , n.ones(Lt/2) )  )

# Triangular:
foo=n.linspace(-1,1,Lt/2,endpoint=False)
Tr_i=n.hstack(  ( foo , foo*-1 )   )

# Dente de Serra:
D_i=n.linspace(-1,1,Lt)

def v(f=220,d=2.,tab=S_i,fv=2.,nu=2.,tabv=S_i):
    if nu==13.789987:
        return n.zeros(int(fa*d))
    Lambda=n.floor(f_a*d)
    ii=n.arange(Lambda)
    Lv=float(len(tabv))

    Gammav_i=n.floor(ii*fv*Lv/f_a) # índices para a LUT
    Gammav_i=n.array(Gammav_i,n.int)
    # padrão de variação do vibrato para cada amostra
    Tv_i=tabv[Gammav_i%int(Lv)] 

    # frequência em Hz em cada amostra
    F_i=f*(   2.**(  Tv_i*nu/12.  )   ) 
    # a movimentação na tabela por amostra
    D_gamma_i=F_i*(Lt/float(f_a))
    Gamma_i=n.cumsum(D_gamma_i) # a movimentação na tabela total
    Gamma_i=n.floor( Gamma_i) # já os índices
    Gamma_i=n.array( Gamma_i, dtype=n.int) # já os índices
    return tab[Gamma_i%int(Lt)] # busca dos índices na tabela

def A(fa=2.,V_dB=10.,d=2.,taba=S_i):
    # Use com: v(d=XXX)*A(d=XXX)
    Lambda=n.floor(f_a*d)
    ii=n.arange(Lambda)
    Lt=float(len(taba))
    Gammaa_i=n.floor(ii*fa*Lt/f_a) # índices para a LUT
    Gammaa_i=n.array(Gammaa_i,n.int)
    # variação da amplitude em cada amostra
    A_i=taba[Gammaa_i%int(Lt)] 
    A_i=1+A_i*(1- 10.**(V_dB/20.))
    return A_i



def adsr(som,A=10.,D=20.,S=-20.,R=100.,xi=1e-2):
    """Envelope ADSR com
    A ataque em milissegundos,
    D decay em milissegundos
    S sustain, com número de decibéis a menos
    R Release em milisegundos

    Atenção para que a duração total é dada pelo som em si
    e que a duração do trecho em sustain é a diferença
    entre a duração total e as durações das partes ADR."""

    a_S=10**(S/20.)
    Lambda=len(som)
    Lambda_A=int(A*f_a*0.001)
    Lambda_D=int(D*f_a*0.001)
    Lambda_R=int(R*f_a*0.001)
    Lambda_S=Lambda - Lambda_A - Lambda_D - Lambda_R

    ii=n.arange(Lambda_A,dtype=n.float)
    A=ii/(Lambda_A-1)
    A_i=A # ok
    ii=n.arange(Lambda_A,Lambda_D+Lambda_A,dtype=n.float)
    D=1-(1-a_S)*(   ( ii-Lambda_A )/( Lambda_D-1) )
    A_i=n.hstack(  (A_i, D  )   )

    S=n.ones(Lambda-Lambda_R-(Lambda_A+Lambda_D),dtype=n.float)*a_S
    A_i=n.hstack( ( A_i, S )  )

    ii=n.arange(Lambda-Lambda_R,Lambda,dtype=n.float)
    R=a_S-a_S*((ii-(Lambda-Lambda_R))/(Lambda_R-1))
    A_i=n.hstack(  (A_i,R)  )

    return som*A_i

triadeM=[0.,4.,7.]
def ac(f=220.,notas=[0.,4.,7.,12.],tab=Q_i,d=2.,nu=0,fv=2.):
    acorde=adsr(v(tab=tab,d=d,f=f*2.**(notas[-1]/12.),nu=nu,fv=fv))
    for na in notas[:-1]:
        acorde+=adsr(v(tab=tab,d=d,f=f*2**(na/12.),nu=nu,fv=fv))
    
    return acorde*10

to= ac(220,[0,4,7])
ms= ac(220,[3,7,10])
ms2=ac(220,[4,8,11])
mi= ac(220,[0,3,8])
mi2=ac(220,[1,4,9])
print("ok 1 - %.2f"%(time.time()-tfoo,)); tfoo=time.time()

to= adsr(ac(220,[0,4,7],nu=2.,fv=10.),R=1500)

ms2=ac(220,[4,8,11])
mi= ac(220,[0,3,8])
mi2=ac(220,[1,4,9])

s=to
s=((s-s.min())/(s.max()-s.min()))*2.-1.
s = n.int16(s * float(2**15-1))
w.write("vomito.wav",f_a, s) # escrita do som

to= adsr(ac(220,[0,4,7],nu=2.,fv=10.),A=60,R=1500)
s=to
s=((s-s.min())/(s.max()-s.min()))*2.-1.
s = n.int16(s * float(2**15-1))
w.write("vomito2.wav",f_a, s) # escrita do som


# Campo harmonico
to= adsr(ac(220,[0,4,7], d=1.,nu=2.,fv=10.), A=260,R=500) 
su= adsr(ac(220,[0,5,9], d=1.,nu=2.,fv=10.), A=260,R=500) 
do= adsr(ac(220,[2,5,11],d=1.,nu=2.,fv=10.),A=260,R=500) 

# medianas
ms= adsr(ac(220,[3,7,10],nu=1.,fv=10.),A=160,R=500)
ms2=adsr(ac(220,[4,8,11],nu=1.,fv=10.),A=160,R=500)
mi=  adsr(ac(220,[0,3,8],nu=1.,fv=10.),A=160,R=500)
mi2= adsr(ac(220,[1,4,9],nu=1.,fv=10.),A=160,R=500)


s=H((to,su,do,to))
s=((s-s.min())/(s.max()-s.min()))*2.-1.
s = n.int16(s * float(2**15-1))
w.write("vomito3.wav",f_a, s) # escrita do som

# Campo harmonico menor
to_m= adsr(ac(220,[0,3,7], d=1.,nu=2.,fv=10.), A=260,R=500) 
su_m= adsr(ac(220,[0,5,8], d=1.,nu=2.,fv=10.), A=260,R=500) 
do_m= adsr(ac(220,[2,5,11],d=1.,nu=2.,fv=10.),A=260,R=500) 

s=H((to_m,su_m,do_m,to_m))
s=((s-s.min())/(s.max()-s.min()))*2.-1.
s = n.int16(s * float(2**15-1))
w.write("vomito4.wav",f_a, s) # escrita do som

# Campo harmonico
to_= adsr(ac(220,[0,4,7], d=1.,nu=0.,fv=10.), A=260,R=500) 
su_= adsr(ac(220,[0,5,9], d=1.,nu=0.,fv=10.), A=260,R=500) 
do_= adsr(ac(220,[2,5,11],d=1.,nu=0.,fv=10.),A=260,R=500) 
do_2= adsr(ac(220,[2,5,7,11],d=1.,nu=0.,fv=10.),A=260,R=500) 

# medianas
ms_ = adsr(ac(220,[3,7,10],nu=0.,fv=10.),A=160,R=500)
ms2_=adsr( ac(220,[4,8,11],nu=0.,fv=10.),A=160,R=500)
mi_ =  adsr(ac(220,[0,3,8],nu=0.,fv=10.),A=160,R=500)
mi2_= adsr( ac(220,[1,4,9],nu=0.,fv=10.),A=160,R=500)


s=H((to_,su_,do_,to_))
s=((s-s.min())/(s.max()-s.min()))*2.-1.
s = n.int16(s * float(2**15-1))
w.write("vomito3b.wav",f_a, s) # escrita do som

s=H((to_,su_,do_,do_2,to_))
s=((s-s.min())/(s.max()-s.min()))*2.-1.
s = n.int16(s * float(2**15-1))
w.write("vomito3b_.wav",f_a, s) # escrita do som


# Campo harmonico menor
to_m_= adsr(ac(220,[0,3,7], d=1.,nu=0.,fv=10.), A=260,R=500) 
su_m_= adsr(ac(220,[0,5,8], d=1.,nu=0.,fv=10.), A=260,R=500) 
do_m_= adsr(ac(220,[2,5,11],d=1.,nu=0.,fv=10.),A=260,R=500) 
do_m_2= adsr(ac(220,[2,5,7,11],d=1.,nu=0.,fv=10.),A=260,R=500) 

s=H((to_m_,su_m_,do_m_,to_m_))
s=((s-s.min())/(s.max()-s.min()))*2.-1.
s = n.int16(s * float(2**15-1))
w.write("vomito4b.wav",f_a, s) # escrita do som

s=H((to_m_,su_m_,do_m_,do_m_2,to_m_))
s=((s-s.min())/(s.max()-s.min()))*2.-1.
s = n.int16(s * float(2**15-1))
w.write("vomito4b_.wav",f_a, s) # escrita do som

s=adsr(v(110.,d=10.,fv=.3)*A(d=10.))
s=((s-s.min())/(s.max()-s.min()))*2.-1.
s = n.int16(s * float(2**15-1))

w.write("chorando.wav",f_a, s) # escrita do som

s=adsr(v(110.,d=10.,fv=.3,tab=D_i)*A(d=10.),S=-5.)
s=((s-s.min())/(s.max()-s.min()))*2.-1.
s = n.int16(s * float(2**15-1))

w.write("chorando2.wav",f_a, s) # escrita do som

s=adsr(v(110.,d=10.,fv=.3,tab=Tr_i)*A(d=10.),S=-5.)
s=((s-s.min())/(s.max()-s.min()))*2.-1.
s = n.int16(s * float(2**15-1))

w.write("chorando3.wav",f_a, s) # escrita do som


s=adsr(v(110.,d=10.,fv=.3,tab=Tr_i)*A(d=10.,fa=.2))
s=((s-s.min())/(s.max()-s.min()))*2.-1.
s = n.int16(s * float(2**15-1))

w.write("chorando4.wav",f_a, s) # escrita do som


s=adsr(v(110.,d=10.,fv=.3,tab=Tr_i)*A(d=10.,fa=12.,V_dB=2.))
s=((s-s.min())/(s.max()-s.min()))*2.-1.
s = n.int16(s * float(2**15-1))

w.write("chorando5.wav",f_a, s) # escrita do som


s=adsr(v(110.,d=10.,fv=.3,tab=Tr_i,tabv=D_i)*A(d=10.,fa=12.,V_dB=2.))
s=((s-s.min())/(s.max()-s.min()))*2.-1.
s = n.int16(s * float(2**15-1))

w.write("chorando6.wav",f_a, s) # escrita do som


s=adsr(v(110.*5,d=10.,fv=.3,tab=Tr_i)*A(d=10.,fa=12.,V_dB=100.),A=100,S=-5.)
s=((s-s.min())/(s.max()-s.min()))*2.-1.
s = n.int16(s * float(2**15-1))

w.write("chorando7.wav",f_a, s) # escrita do som




s=adsr(v(110.*5,d=10.,fv=.3,tab=Tr_i)*A(d=10.,fa=12.,V_dB=100.),A=100,S=-5.)
s=((s-s.min())/(s.max()-s.min()))*2.-1.
s = n.int16(s * float(2**15-1))

w.write("chorando7.wav",f_a, s) # escrita do som

####
# ruidos
Lambda = 100000  # Lambda sempre par
# diferença das frequências entre coeficiêntes vizinhos:
df = f_a/float(Lambda)


coefs = n.exp(1j*n.random.uniform(0, 2*n.pi, Lambda))
# real par, imaginaria impar
coefs[Lambda/2+1:] = n.real(coefs[1:Lambda/2])[::-1] - 1j * \
    n.imag(coefs[1:Lambda/2])[::-1]
coefs[0] = 0.  # sem bias
coefs[Lambda/2] = 1.  # freq max eh real simplesmente

# as frequências relativas a cada coeficiente
# acima de Lambda/2 nao vale
fi = n.arange(coefs.shape[0])*df
f0 = 15.  # iniciamos o ruido em 15 Hz
i0 = n.floor(f0/df)  # primeiro coef a valer
coefs[:i0] = n.zeros(i0)
f0 = fi[i0]

# obtenção do ruído em suas amostras temporais
ruido = n.fft.ifft(coefs)
r = n.real(ruido)
r = ((r-r.min())/(r.max()-r.min()))*2-1
r_=r
r = n.int16(r * float(2**15-1))
w.write('branco.wav', f_a, r)

def N(arr):
    r=arr
    r = ((r-r.min())/(r.max()-r.min()))*2-1
    return n.int16(r * float(2**15-1))
    

fator = 10.**(-6/20.)
alphai = fator**(n.log2(fi[i0:]/f0))
c = n.copy(coefs)
c[i0:] = c[i0:]*alphai

# real par, imaginaria impar
c[Lambda/2+1:] = n.real(c[1:Lambda/2])[::-1] - 1j * \
    n.imag(c[1:Lambda/2])[::-1]

# realizando amostras temporais do ruído marrom
ruido = n.fft.ifft(c)
r = n.real(ruido)
r = ((r-r.min())/(r.max()-r.min()))*2-1
r__=r
r = n.int16(r * float(2**15-1))
w.write('marrom.wav', f_a, r)

### 2.53 Ruído azul
# para cada oitava, ganhamos 3dB
fator = 10.**(3/20.)
alphai = fator**(n.log2(fi[i0:]/f0))
c = n.copy(coefs)
c[i0:] = c[i0:]*alphai

# real par, imaginaria impar
c[Lambda/2+1:] = n.real(c[1:Lambda/2])[::-1] - 1j * \
    n.imag(c[1:Lambda/2])[::-1]

# realizando amostras temporais do ruído azul
ruido = n.fft.ifft(c)
r = n.real(ruido)
ra=r
r = ((r-r.min())/(r.max()-r.min()))*2-1
r = n.int16(r * float(2**15-1))
w.write('azul.wav', f_a, r)

### 2.54 Ruido violeta
# a cada oitava, ganhamos 6dB
fator = 10.**(6/20.)
alphai = fator**(n.log2(fi[i0:]/f0))
c = n.copy(coefs)
c[i0:] = c[i0:]*alphai

# real par, imaginaria impar
c[Lambda/2+1:] = n.real(c[1:Lambda/2])[::-1] - 1j * \
    n.imag(c[1:Lambda/2])[::-1]

ruido = n.fft.ifft(c)
r = n.real(ruido)
rv=r
r = ((r-r.min())/(r.max()-r.min()))*2-1
r = n.int16(r * float(2**15-1))
w.write('violeta.wav', f_a, r)




w.write('fedendo.wav', f_a, H((r,r,r,r,r)))
w.write('fedendo2.wav', f_a, N(H((r_,r_,r_,r_,r_))))


w.write('varrendo.wav', f_a, N(H((adsr(r_[:int(.5*f_a)],A=250),n.zeros(.3*f_a),
                                  adsr(r_[:int(.5*f_a)],A=250),n.zeros(.3*f_a),
                                  adsr(r_[:int(.5*f_a)],A=250),n.zeros(.3*f_a),
                                  adsr(r_[:int(.5*f_a)],A=250),n.zeros(.3*f_a),
                                  adsr(r_[:int(.5*f_a)],A=250),n.zeros(.3*f_a),
                                  adsr(r_[:int(.5*f_a)],A=250)))))

w.write('varrendo2.wav', f_a, N(H((adsr(r__[:int(.5*f_a)],A=250),n.zeros(.3*f_a),
                                   adsr(r__[:int(.5*f_a)],A=250),n.zeros(.3*f_a),
                                   adsr(r__[:int(.5*f_a)],A=250),n.zeros(.3*f_a),
                                   adsr(r__[:int(.5*f_a)],A=250),n.zeros(.3*f_a),
                                   adsr(r__[:int(.5*f_a)],A=250),n.zeros(.3*f_a),
                                   adsr(r__[:int(.5*f_a)],A=250)))))


w.write('varrendo3.wav', f_a, N(H((adsr(ra[:int(.5*f_a)],A=250),n.zeros(.3*f_a),
                                   adsr(ra[:int(.5*f_a)],A=250),n.zeros(.3*f_a),
                                   adsr(ra[:int(.5*f_a)],A=250),n.zeros(.3*f_a),
                                   adsr(ra[:int(.5*f_a)],A=250),n.zeros(.3*f_a),
                                   adsr(ra[:int(.5*f_a)],A=250),n.zeros(.3*f_a),
                                   adsr(ra[:int(.5*f_a)],A=250)))))

varr=N(H((adsr(rv[:int(.5*f_a)],A=250),n.zeros(.3*f_a),
                                   adsr(rv[:int(.5*f_a)],A=250),n.zeros(.3*f_a),
                                   adsr(rv[:int(.5*f_a)],A=250),n.zeros(.3*f_a),
                                   adsr(rv[:int(.5*f_a)],A=250),n.zeros(.3*f_a),
                                   adsr(rv[:int(.5*f_a)],A=250),n.zeros(.3*f_a),
                                   adsr(rv[:int(.5*f_a)],A=250))))

w.write('varrendo4.wav', f_a, varr)
w.write('varrendo5.wav', f_a, varr[::3])
w.write('varrendo6.wav', f_a, varr[::4])
w.write('cocando.wav',   f_a, varr[::5])
w.write('cocando2.wav',  f_a, varr[::6])
w.write('cocando3.wav',  f_a, varr[::7])
w.write('cocando4.wav',  f_a, varr[::8])
w.write('cocando5.wav',  f_a, N(rv[:2*f_a]*A(d=2.)))
w.write('cocando6.wav',  f_a, N(rv[:2*f_a]*A(d=2.,taba=D_i[::-1],fa=5.)))
w.write('cocando7.wav',  f_a, N(rv[:2*f_a]*A(d=2.,taba=Q_i[::-1],fa=5.)))
w.write('cocando8.wav',  f_a, N(ra[:2*f_a]*A(d=2.,taba=D_i[::-1],fa=15.)))
w.write('cocando9.wav',  f_a, N(rv[:2*f_a]*A(d=2.,taba=D_i[::-1],fa=15.)))

w.write('pulando.wav',  f_a, N(H((adsr(v(d=.5,tab=Q_i)),
                               adsr(v(d=.5,tab=Q_i)),
                               adsr(v(d=.5,tab=Q_i)),
                               adsr(v(d=.5,tab=Q_i)),
                                ))))

w.write('pulando2.wav',  f_a, N(H((adsr(v(d=.5,tab=Q_i,tabv=D_i)),
                                   adsr(v(d=.5,tab=Q_i,tabv=D_i)),
                                   adsr(v(d=.5,tab=Q_i,tabv=D_i)),
                                   adsr(v(d=.5,tab=Q_i,tabv=D_i)),
                                ))))



w.write('pulando3.wav',  f_a, N(H((adsr(v(d=.5,tab=Q_i,nu=24.,tabv=D_i)),
                                   adsr(v(d=.5,tab=Q_i,nu=24.,tabv=D_i)),
                                   adsr(v(d=.5,tab=Q_i,nu=24.,tabv=D_i)),
                                   adsr(v(d=.5,tab=Q_i,nu=24.,tabv=D_i)),
                                ))))


w.write('pulando4.wav',  f_a, N(H((adsr(v(f=880,d=.5,tab=Q_i,nu=24.,tabv=D_i)),
                                   adsr(v(f=880,d=.5,tab=Q_i,nu=24.,tabv=D_i)),
                                   adsr(v(f=880,d=.5,tab=Q_i,nu=24.,tabv=D_i)),
                                   adsr(v(f=880,d=.5,tab=Q_i,nu=24.,tabv=D_i)),
                                ))))


w.write('pulando5.wav',  f_a, N(H((adsr(v(f=880,d=.5,tab=D_i,nu=4.,tabv=D_i)),
                                   adsr(v(f=880,d=.5,tab=D_i,nu=4.,tabv=D_i)),
                                   adsr(v(f=880,d=.5,tab=D_i,nu=4.,tabv=D_i)),
                                   adsr(v(f=880,d=.5,tab=D_i,nu=4.,tabv=D_i)),
                                ))))

w.write('masturba.wav',  f_a, N(H((
                                   adsr(r_[:int(.100*f_a)],R=30),n.zeros(int(f_a*.1)),
                                   adsr(r_[:int(.100*f_a)],R=30),n.zeros(int(f_a*.1)),
                                   adsr(r_[:int(.100*f_a)],R=30),n.zeros(int(f_a*.1)),
                                   adsr(r_[:int(.100*f_a)],R=30),n.zeros(int(f_a*.1)),
                                ))))

w.write('masturba2.wav',  f_a, N(H((
                                   adsr(r_[:int(.200*f_a)]),n.zeros(int(f_a*.1)),
                                   adsr(r_[:int(.200*f_a)]),n.zeros(int(f_a*.1)),
                                   adsr(r_[:int(.200*f_a)]),n.zeros(int(f_a*.1)),
                                   adsr(r_[:int(.200*f_a)]),n.zeros(int(f_a*.1)),
                                ))))

w.write('masturba3.wav',  f_a, N(H((
                                   adsr(r_[:int(.02*f_a)],A=5,D=5,R=10),n.zeros(int(f_a*.1)),
                                   adsr(r_[:int(.02*f_a)],A=5,D=5,R=10),n.zeros(int(f_a*.1)),
                                   adsr(r_[:int(.02*f_a)],A=5,D=5,R=10),n.zeros(int(f_a*.1)),
                                   adsr(r_[:int(.02*f_a)],A=5,D=5,R=10),n.zeros(int(f_a*.1)),
                                ))))

w.write('masturba4.wav',  f_a,  N(adsr(r_, A=5.,D=5.,R=10.)*A(d=len(r_)/f_a,fa=15)))
w.write('masturba4b.wav',  f_a, N(adsr(r__,A=5.,D=5.,R=10.)*A(d=len(r_)/f_a,fa=15)))
w.write('masturba4c.wav',  f_a, N(adsr(ra ,A=5.,D=5.,R=10.)*A(d=len(r_)/f_a,fa=15)))
w.write('masturba4d.wav',  f_a, N(adsr(rv ,A=5.,D=5.,R=10.)*A(d=len(r_)/f_a,fa=15)))

w.write('masturba5.wav',  f_a,  N(adsr(r_, A=5.,D=5.,R=10.)*A(d=len(r_)/f_a,fa=25)))
w.write('masturba5b.wav',  f_a, N(adsr(r__,A=5.,D=5.,R=10.)*A(d=len(r_)/f_a,fa=25)))
w.write('masturba5c.wav',  f_a, N(adsr(ra ,A=5.,D=5.,R=10.)*A(d=len(r_)/f_a,fa=25)))
w.write('masturba5d.wav',  f_a, N(adsr(rv ,A=5.,D=5.,R=10.)*A(d=len(r_)/f_a,fa=25)))

w.write('masturba6.wav',  f_a,  N(adsr(r_, A=5.,D=5.,R=10.)*A(d=len(r_)/f_a,fa=25,V_dB=50.)))
w.write('masturba6b.wav',  f_a, N(adsr(r__,A=5.,D=5.,R=10.)*A(d=len(r_)/f_a,fa=25,V_dB=50.)))
w.write('masturba6c.wav',  f_a, N(adsr(ra ,A=5.,D=5.,R=10.)*A(d=len(r_)/f_a,fa=25,V_dB=50.)))
w.write('masturba6d.wav',  f_a, N(adsr(rv ,A=5.,D=5.,R=10.)*A(d=len(r_)/f_a,fa=25,V_dB=50.)))


w.write('masturba7.wav',  f_a,  N(adsr(r_, A=5.,D=5.,R=10.,S=-5.)*A(d=len(r_)/f_a,fa=5,V_dB=50.)))
w.write('masturba7b.wav',  f_a, N(adsr(r__,A=5.,D=5.,R=10.,S=-5.)*A(d=len(r_)/f_a,fa=5,V_dB=50.)))
w.write('masturba7c.wav',  f_a, N(adsr(ra ,A=5.,D=5.,R=10.,S=-5.)*A(d=len(r_)/f_a,fa=5,V_dB=50.)))
w.write('masturba7d.wav',  f_a, N(adsr(rv ,A=5.,D=5.,R=10.,S=-5.)*A(d=len(r_)/f_a,fa=5,V_dB=50.)))






sys.exit()
#frase="semicondutor livre"
#arq=frase.split()[0]
#os.system("espeak -vpt-pt+%s -w%s.wav '%s'"%(random.sample(vozes,1)[0],arq,frase))
#ff=w.read("%s.wav"%(arq,))[1]
#ff_=n.fft.fft(ff)
#s=ff2=n.fft.ifft( n.hstack((ff_,n.zeros(len(ff_)) ))  ).real
#sc_aud=((s-s.min())/(s.max()-s.min()))*2.-1.

vozes="f3,f2,f1,f5,m5,m1,m3".split(",")
def fala(frase="Semicondutor livre"):
    arq=frase.split()[0]
    #os.system("espeak -vpt-pt+%s -w%s.wav -g110 -p99 -s110 -b=1 '%s'"%(random.sample(vozes,1)[0],arq,frase))
    os.system("espeak -vpt-pt+%s -w%s.wav -p99 -b=1 '%s'"%(random.sample(vozes,1)[0],arq,frase))
    #os.system("espeak -vpt-pt+%s -w%s.wav -g110 -p99 -s130 -b=1 '%s'"%(random.sample(vozes,1)[0],arq,frase))
    ff=w.read("%s.wav"%(arq,))[1]
    ff_=n.fft.fft(ff)
    s=ff2=n.fft.ifft( n.hstack((ff_,n.zeros(len(ff_)) ))  ).real
    sc_aud=((s-s.min())/(s.max()-s.min()))*2.-1.
    return sc_aud*10

mano=fala("maaaannnnnnoooooooo")
mano2=fala("mannnnnooooooooooooooooooooooooooooo")
lm=fala(u"eh o labi macambira")
T=f_a # 60 BPM

#linha=H((mano,n.zeros(f_a-len(mano)),)
linha=n.zeros((4*f_a)*16) # 16 compassos 4/4 a 60 BPM

# COMPASSO 1
linha[:len(mano)]=mano
TT=2*f_a; som=mano2
linha[TT:TT+len(som)]+=som

# COMPASSO 2
# o segundo compasso é silêncio que fecha com a fala lm
linha[8*f_a-len(lm):8*f_a]+=lm

# COMPASSO 3
fb=110
pb=[(0,1),(7,1),(0,1),(7,.75),(7,.25)]
lb=H([v(f=fb*2**(i[0]/12.),d=i[1],nu=0,tab=D_i) for i in pb])

fm=300
pm=[(0,.5),(-1,.5),(4,.5),(-7,.5),
    (0,.5),(11,.5),(12,.5),(2,.5)]
lm=H([v(f=fm*2**(i[0]/12.),d=i[1],fv=15.,nu=.5,tab=Tr_i) for i in pm])

fa=2200
pa=[(0,3),(2,1)]
la=H([v(f=fa*2**(i[0]/12.),d=i[1],nu=2,fv=1,tab=S_i) for i in pa])

TT=8*f_a; som=lb
linha[TT:TT+len(som)]+=som
TT=8*f_a; som=lm
linha[TT:TT+len(som)]+=som
TT=8*f_a; som=la
linha[TT:TT+len(som)]+=som
#linha[TT:TT+len(som)]*=(1/3.)

# COMPASSO 4
pa=[(4,2),(5,2)]
la=H([v(f=fa*2**(i[0]/12.),d=i[1],nu=3,fv=1,tab=S_i,tabv=D_i) for i in pa])

TT=12*f_a; som=lb
linha[TT:TT+len(som)]+=som
TT=12*f_a; som=lm
linha[TT:TT+len(som)]+=som
TT=12*f_a; som=la
linha[TT:TT+len(som)]+=som
linha[TT:TT+len(som)]*=(1/3.)

# COMPASSO 5: acabou a apresentação do material,
# inicia desenvolvimento,
# só a linha do grave comeca a fazer firula
#fb=110
pb=[(0,1,4.,1.),(7,1,4.,2.),
    (0,1,4.,2.),(7,.75,8.,2.),(7,.25,8.,16.)]
lb=H([v(f=fb*2**(i[0]/12.),d=i[1],fv=i[2],nu=i[3],tab=D_i) for i in pb])
TT=16*f_a; som=lb
linha[TT:TT+len(som)]+=som

# COMPASSO 6: linha grave continua firulando
# agora com uma continuacao da linha
# e fala lm no segundo tempo
pb=[(0,1,1.,24.),(-1,1,4.,2.),
    (-7,1,14.,2.),(7,.75,8.,2.),(11,.25,108.,16.)]
lb=H([v(f=fb*2**(i[0]/12.),d=i[1],fv=i[2],nu=i[3],tab=Q_i) for i in pb])
TT=20*f_a; som=lb
linha[TT:TT+len(som)]+=som
TT=21*f_a; som=fala("eh o laaaaaab maaaaaacambira a a a aa a")
linha[TT:TT+len(som)]+=som

# COMPASSO 7: as 3 linhas voltam, fazem certo furduncio
pb=[(0,1,4.,1.),(7,1,4.,2.),
    (0,1,104.,2.),(7,.75,8.,2.),(7,.25,8.,16.)]
lb=H([v(f=fb*2**(i[0]/12.),d=i[1],fv=i[2],nu=i[3],tab=D_i) for i in pb])

fm=300
pm=[(0,.5,3.,3.),(-1,.5,0,0),(4,.5,2.,2.),(-7,.5,0,0),
    (0,.5,1.,1.),(11,.5,0,0),(12,.5,104.,2.),(2,.5,0,0)]
lm=H([v(f=fm*2**(i[0]/12.),d=i[1],fv=15.,nu=.5,tab=Tr_i) for i in pm])

print("ok 2 - %.2f"%(time.time()-tfoo,)); tfoo=time.time()
fa=2200
pa=[(0,3,1,4),(2,1,10,4)]
ab=H([v(f=fa*2**(i[0]/12.),d=i[1],nu=i[3],fv=i[2],tab=S_i) for i in pa])

TT=6*(4*f_a); som=lb # passaram-se seis compassos
linha[TT:TT+len(som)]+=som
TT=6*(4*f_a); som=lm
linha[TT:TT+len(som)]+=som
TT=6*(4*f_a); som=la
linha[TT:TT+len(som)]+=som

# COMPASSO 8: dev

pb=[(0,1,4.,1.),(7,1,4.,2.),
    (0,1,104.,2.),(7,.75,8.,2.),(7,.25,8.,16.)]
lb=H([v(f=fb*2**(i[0]/12.),d=i[1],fv=i[2],nu=i[3],tab=S_i) for i in pb])

fm=300
pm=[(0,.5,3.,3.),(-1,.5,0,0),(4,.5,2.,2.),(-7,.5,0,0),
    (0,.5,1.,1.),(11,.5,0,0),(12,.5,104.,2.),(2,.5,0,0)]
lm=H([v(f=fm*2**(i[0]/12.),d=i[1],fv=15.,nu=.5,tab=S_i) for i in pm])

fa=2200
pa=[(0,1,110,4),(2,1,110,40),(0,1,10,40),(2,1,10,4)]
ab=H([v(f=fa*2**(i[0]/12.),d=i[1],nu=i[3],fv=i[2],tab=S_i) for i in pa])

TT=7*(4*f_a); som=lb # passaram-se seis compassos
linha[TT:TT+len(som)]+=som
TT=7*(4*f_a); som=lm
linha[TT:TT+len(som)]+=som
TT=7*(4*f_a); som=la
linha[TT:TT+len(som)]+=som

# COMPASSO 9: encadeamentos harmonicos

to=ac(110,[0,4,7,12],tab=S_i)
sr=ac(110,[2,5,9,14],tab=S_i)

TT=8*(4*f_a); som=to
linha[TT:TT+len(som)]+=som
TT=8*(4*f_a)+2*f_a; som=sr
linha[TT:TT+len(som)]+=som
# e das linhas, entrecortadas

pb=[(0,1,0.,0.),(7,1,0.,0.),
    (0,1,0.,13.789987),(7,.75,0.,13.789987),(7,.25,0.,0.)]
lb=H([v(f=fb*2**(i[0]/12.),d=i[1],fv=i[2],nu=i[3],tab=S_i) for i in pb])

fm=300
pm=[(0,.5,3.,13.789987),(-1,.5,0,0),(4,.5,2.,13.789987),(-7,.5,0,0),
    (0,.5,1.,13.789987),(11,.5,0,13.789987),(12,.5,104.,2.),(2,.5,500.,0)]
lm=H([v(f=fm*2**(i[0]/12.),d=i[1],fv=15.,nu=.5,tab=S_i) for i in pm])

fa=2200
pa=[(0,1,110,4),(2,1,110,13.789987),(0,1,10,4),(2,1,10,14)]
la=H([v(f=fa*2**(i[0]/12.),d=i[1],nu=i[3],fv=i[2],tab=S_i) for i in pa])

TT=8*(4*f_a); som=lb # passaram-se oito compassos
linha[TT:TT+len(som)]+=som
TT=8*(4*f_a); som=lm
linha[TT:TT+len(som)]+=som
TT=8*(4*f_a); som=la
linha[TT:TT+len(som)]+=som

TT=8*(4*f_a); som=fala("mu u u u ui i i itoo o o o ma a a assa maa a a aa a nnnnnooooo")
linha[TT:TT+len(som)]+=som
# COMPASSO 10, dominante
d7=ac(110,[2,5,7,11],tab=S_i)

d7_=ac(110,[2,5,7,11,-1,-5,14,19],tab=S_i)
TT=9*(4*f_a); som=d7
linha[TT:TT+len(som)]+=som
TT=9*(4*f_a)+2*f_a; som=d7_
linha[TT:TT+len(som)]+=som


pb=[(0,1,0.,13.789987),(7,1,0.,13.789987),
    (0,1,0.,0),(7,.75,0.,0),(7,.25,0.,13.789987)]
lb=H([v(f=fb*2**(i[0]/12.),d=i[1],fv=i[2],nu=i[3],tab=S_i) for i in pb])

fm=300
pm=[(0,.5,3.,0.),(-1,.5,0,13.789987),(4,.5,2.,0),(-7,.5,0,13.789987),
    (0,.5,1.,0),(11,.5,0,0),(12,.5,104.,13.789987),(2,.5,500.,13.789987)]
lm=H([v(f=fm*2**(i[0]/12.),d=i[1],fv=15.,nu=.5,tab=S_i) for i in pm])

fa=2200
pa=[(0,1,110,13.789987),(2,1,110,0),(0,1,10,400),(2,1,10,14)]
la=H([v(f=fa*2**(i[0]/12.),d=i[1],nu=i[3],fv=i[2],tab=S_i) for i in pa])

TT=9*(4*f_a); som=lb # passaram-se oito compassos
linha[TT:TT+len(som)]+=som
TT=9*(4*f_a); som=lm
linha[TT:TT+len(som)]+=som
TT=9*(4*f_a); som=la
linha[TT:TT+len(som)]+=som
TT=int(9.5*(4*f_a)); som=fala("mu u u u ui i i itoo o o o ma a a assa maa a a aa a nnnnnooooo")
linha[TT:TT+len(som)]+=som

# COMPASSO 11, volta p tônica e modula para a dominante

to=ac(110,[0,4,7,12],tab=Tr_i,d=1.)
sr=ac(110,[2,5,9,14],tab=Tr_i,nu=.2,d=1.)
do=ac(110,[2,-1,7,11],tab=Tr_i,d=1.)
d_do=ac(110,[2,6,9,12],tab=S_i,nu=.2,d=1.)
print("ok 3 - %.2f"%(time.time()-tfoo,)); tfoo=time.time()

TT=10*(4*f_a); som=to
linha[TT:TT+len(som)]+=som
TT=10*(4*f_a)+f_a; som=sr
linha[TT:TT+len(som)]+=som
TT=10*(4*f_a)+2*f_a; som=do
linha[TT:TT+len(som)]+=som
TT=10*(4*f_a)+3*f_a; som=d_do
linha[TT:TT+len(som)]+=som

fa=2200
pa=[(0.,1.,1.,.2),(2.,1.,2.,.2),(4.,1.,1.,.2),(5.,1.,2.,14.)]
la=H([v(f=fa*2**(i[0]/12.),d=i[1],nu=i[3],fv=i[2],tab=S_i) for i in pa])

TT=10*(4*f_a); som=la
linha[TT:TT+len(som)]+=som
#TT=int(10.75*(4*f_a)); som=fala("python, javascripiti, peagahpe, dijango, bashhhhhshshshs, linuquis, a agah teee, aaaa aaaagggaaaa, teeee, freaaaaqueee, coooddddiiiinnnnguee, ehhheeeehhehe nnoiiiiissss maccaaaammbbiiiiiiiirrrraaaa  aaa  a a a a  a a a aaaaaaaaaaa  aaaaaaaa raaaaaaaataaaaaaaa taaaaaa raaaa tataata raraaraiaiaiaiaoaoaooaapapapaparararrawawawaagagaggafafafajajajaoaoabababadadacacaccaa")
TT=int(10.75*(4*f_a)); som=fala("python, javascripiti, peagahpe, dijango, bashhhhhshshshs, linuquis, a agah teee, aaaa aaaagggaaaa, teeee, freaaaaqueee, coooddddiiiinnnnguee")
linha[TT:TT+len(som)]+=som

# COMPASSO 12: ápice, dominante e seu s espaço harmonico
do=ac(110,[2,-1,7,11],tab=S_i,d=1.)
s_do=ac(110,[-3,0,4,9],tab=S_i,d=1.)
d_do=ac(110,[2,6,9,12],tab=S_i,nu=.2)

TT=11*(4*f_a); som=do
linha[TT:TT+len(som)]+=som
TT=11*(4*f_a)+f_a; som=s_do
linha[TT:TT+len(som)]+=som
TT=11*(4*f_a)+2*f_a; som=d_do
linha[TT:TT+len(som)]+=som

fa=2200
pa=[(7.,1.,1.,.2),(9.,1.,2.,.2),(7.,1.,1.,.2),(9.,1.,2.,14.)]
la=H([v(f=fa*2**(i[0]/12.),d=i[1],nu=i[3],fv=i[2],tab=S_i) for i in pa])

fm=300
pm=[(7,.5,3.,0.),(6,.5,0,13.789987),(11,.5,2.,0),(2,.5,0,13.789987),
    (7,.5,1.,0),(6,.5,0,0),(7,.5,104.,13.789987),(9,.5,500.,5.)]
lm=H([v(f=fm*2**(i[0]/12.),d=i[1],fv=15.,nu=.5,tab=S_i) for i in pm])
TT=11*(4*f_a); som=lm
linha[TT:TT+len(som)]+=som
TT=11*(4*f_a); som=lm
linha[TT:TT+len(som)]+=som

# COMPASSO 13
# dominante, dominante da dominante com setima
# baixa meio tom, vira d7 sub
do=ac(110,[2,-1,7,11],tab=S_i,d=1.)
s_do=ac(110,[-3,0,4,9],tab=S_i,d=1.)
d_do=ac(110,[2,6,9,12],tab=S_i,nu=.1)
d_sub=ac(110,[1,5,8,11],tab=S_i,nu=.1)

TT=12*(4*f_a); som=do
linha[TT:TT+len(som)]+=som
TT=12*(4*f_a)+f_a; som=s_do
linha[TT:TT+len(som)]+=som
TT=12*(4*f_a)+2*f_a; som=d_do
linha[TT:TT+len(som)]+=som
TT=12*(4*f_a)+3*f_a; som=d_sub
linha[TT:TT+len(som)]+=som

# COMPASSO 14, volta p tônica, sossega

do=ac(110,[2,-1,7,11],  tab=Tr_i,d=1.)
s_do=ac(110,[-3,0,4,9], tab=Tr_i,d=1.)
d_do=ac(110,[2,6,9,12], tab=Tr_i,nu=.1)
d_sub=ac(110,[1,5,8,11],tab=Tr_i,nu=.1)

TT=13*(4*f_a); som=do
linha[TT:TT+len(som)]+=som
TT=13*(4*f_a)+f_a; som=s_do
linha[TT:TT+len(som)]+=som
TT=13*(4*f_a)+2*f_a; som=d_do
linha[TT:TT+len(som)]+=som
TT=13*(4*f_a)+3*f_a; som=d_sub
linha[TT:TT+len(som)]+=som

fm=13000
pm=[(7,.5,3.,0.),(-6,.5,0,13.789987),(-11,.5,2.,0),(-2,.5,0,13.789987),
    (-27,.5,1.,0),(-36,.5,0,0),(-47,.5,104.,2),(-59,.5,500.,5.)]
lm=H([v(f=fm*2**(i[0]/12.),d=i[1],fv=15.,nu=.5,tab=S_i) for i in pm])
TT=13*(4*f_a); som=lm
linha[TT:TT+len(som)]+=som

# COMPASSO 15
# faz cadencia I-III-V-I p firmar tonica

to=ac(110,[0,4,7,12],tab=Tr_i,d=1.)
sr=ac(110,[2,5,9,14],tab=Tr_i,nu=.2,d=1.)
d_sub=ac(110,[1,5,8,11],tab=Tr_i,nu=.1)

TT=14*(4*f_a); som=do
linha[TT:TT+len(som)]+=som
TT=14*(4*f_a)+f_a; som=sr
linha[TT:TT+len(som)]+=som
TT=14*(4*f_a)+2*f_a; som=d_sub
linha[TT:TT+len(som)]+=som

print("ok 4 - %.2f"%(time.time()-tfoo,)); tfoo=time.time()
TT=int(14.25*(4*f_a)); som=fala("democracia direta, economia solidaria, autogestao coletiva, agora comuououououns, irmandade do codigooooo")
linha[TT:TT+len(som)]+=som
# COMPASSO 16: termina
# firula no acorde de tonica
# com falas ou outras coisas

t= v(f=fb*2**(0/12.),  d=1.,nu=0,fv=0,tab=S_i)
tr=v(f=fb*2**(4/12.), d=1.,nu=0,fv=0,tab=S_i)
d= v(f=fb*2**(7/12.),  d=1.,nu=0,fv=0,tab=S_i)
t_=v(f=fb*2**(12/12.),d=1.,nu=0,fv=0,tab=S_i)


t2= v(f=fb*2**(0/12.) ,d=1.,nu=2.,fv=2.,tab=S_i)
tr2=v(f=fb*2**(4/12.),d=1.,nu=2.,fv=2.,tab=S_i)
d2= v(f=fb*2**(7/12.) ,d=1.,nu=4,fv=4.,tab=S_i)

TT=15*(4*f_a)

TT=15*(4*f_a)+0*f_a; som=t
linha[TT:TT+len(som)]+=som
TT=15*(4*f_a)+0*f_a; som=d
linha[TT:TT+len(som)]+=som

TT=15*(4*f_a)+1*f_a; som=t
linha[TT:TT+len(som)]+=som
TT=15*(4*f_a)+1*f_a; som=d2
linha[TT:TT+len(som)]+=som
TT=15*(4*f_a)+1*f_a; som=tr
linha[TT:TT+len(som)]+=som

TT=15*(4*f_a)+2*f_a; som=t
linha[TT:TT+len(som)]+=som
TT=15*(4*f_a)+2*f_a; som=d2
linha[TT:TT+len(som)]+=som
TT=15*(4*f_a)+2*f_a; som=tr
linha[TT:TT+len(som)]+=som
TT=15*(4*f_a)+2*f_a; som=t2
linha[TT:TT+len(som)]+=som
TT=15*(4*f_a)+2*f_a; som=t_
linha[TT:TT+len(som)]+=som

TT=15*(4*f_a)+3*f_a; som=t
linha[TT:TT+len(som)]+=som
TT=15*(4*f_a)+3*f_a; som=tr2
linha[TT:TT+len(som)]+=som
TT=15*(4*f_a)+3*f_a; som=t_
linha[TT:TT+len(som)]+=som
linha=H((linha,v(f=fb*2**(0/12.),  d=4.,nu=0,fv=0,tab=S_i)+v(f=fb*2**(16/12.),  d=4.,nu=.2,fv=1,tab=S_i)+v(f=fb*2**(24/12.),  d=4.,nu=1.5,fv=6,tab=S_i)))

print("ok 5 - %.2f"%(time.time()-tfoo,)); tfoo=time.time()
s=linha
s=((s-s.min())/(s.max()-s.min()))*2.-1.
s = n.int16(s * float(2**15-1))
w.write("penalva.wav",f_a, s) # escrita do som



sys.exit()

f1=fala("semicondutor livre")
TT=f_a*2*2
s[TT:TT+len(f1)]=f1*2
print("ok 2")
f2=fala(u"Cabreuhbo Peixoto")
TT=f_a*5*2
s[TT:TT+len(f2)]=f2*2
f3=fala(u"Fim do mundo, eh o Maia descendente do Nostradamus")
TT=f_a*8*2
s[TT:TT+len(f3)]=f3*2
f4=fala(u"Comofaz no espasso sideral")
print("ok 3")
TT=f_a*11*2
s[TT:TT+len(f4)]=f4*2
print("ok 3b")
f5=fala(u"Poh de estrela, eh a Yupana construindo castelo de areia")
print("ok 4a")
TT=int(f_a*11.5*2)
s[TT:TT+len(f5)]=f5*2
print("ok 4")

f5=fala(u"abaco")
TT=int(f_a*17.5*2)
s[TT:TT+len(f5)]=f5*2
f5=fala(u"oraculo")
TT=int(f_a*21.5*2)
s[TT:TT+len(f5)]=f5*2

print("ok 5")

s=((s-s.min())/(s.max()-s.min()))*2.-1.
s = n.int16(s * float(2**15-1))
w.write("soares.wav",f_a, s) # escrita do som
sys.exit()

s=ac(220.,triadeM)
s2=ac(220.,triadeM,Tr_i)
#s=n.hstack(( s,intervaloHarmonico(220,i) ))
s=H((s,s2,s,s2))


s1=ac(220.,triadeM)
s2=ac(220.,[0,5,9]) # subdominante
s3=ac(220.,[2,7,11]) # dominante
s4=ac(220.,[2,5,9]) # sub relativa
s5=ac(220.,[0,4,9]) # ton relat / sub anti


s=H((s,s1,s2,s3,s1,  s1,s4,s3,s1,  s1,s2,s3,s5, s5,s4,s3,s1  ))
s_=n.copy(s)


s=((s-s.min())/(s.max()-s.min()))*2.-1.

# most music players read only 16-bit wav files, so let's convert the array
s = n.int16(s * float(2**15-1))
w.write("acordeCedo.wav",f_a, s) # escrita do som

f0=110.
padrao=[0,0,-1,0]
padrao2=[10,0,-1,50]
padrao3=[10,0,-1,50]
padrao4=[10,0,-5,5]
padrao5=[11,0,-1,7]
pf=[f0*2**(pi/12.) for pi in padrao]
som1=n.hstack([adsr(v(f=ff,d=0.5,tab=Q_i)) for ff in pf])
som2=n.hstack([adsr(v(f=f0*2**(pi/12.),d=0.5,tab=Q_i)) for pi in padrao2])
som3=n.hstack([adsr(v(f=f0*2**(pi/12.),d=0.5,tab=Q_i)) for pi in padrao3])
som4=n.hstack([adsr(v(f=f0*2**(pi/12.),d=0.5,tab=Q_i)) for pi in padrao4])
som5=n.hstack([adsr(v(f=f0*2**(pi/12.),d=0.5,tab=Q_i)) for pi in padrao5])

s=n.hstack((som1,som2,som5,som1,som4,som2,som3,som1))
s__=n.copy(s)
s=((s-s.min())/(s.max()-s.min()))*2.-1.
s = n.int16(s * float(2**15-1))
w.write("linha.wav",f_a, s) # escrita do som


# para sobrepor com s_
# Intro:
padrao=[0,0,-1,0]
padrao1=[0,4,7,12]
padrao2=[0,-5,-9,-12]
padrao3=[-9,-5,-3,-1]
pp=padrao+padrao1+padrao2+padrao3
som=n.hstack([adsr(v(f=f0*2**(pi/12.),d=0.5,tab=Q_i)) for pi in pp])

#Ciclo1
padrao1=[(i,.5) for i in [0,12,24,36]]
padrao2=[(0,.25),(0,.25),(0,.5),(24+5,.5),(36+9,.5)]
padrao3=[(0,.25),(-1,.25),(-5,.5),(24+7,.5),(36+11,.5)]
padrao4_1=[(0,.25),(0,.25),(0,.5),(-24,.5),(36*2,.5)]
pp_=padrao1+padrao2+padrao3+padrao4_1

#Ciclo2
padrao1=[(i,.5) for i in [0,12,24,36]]
padrao2=[(0,.25),(2,.25),(2,.5),(24+2,.5),(36+9,.5)]
padrao3=[(0,.25),(-1,.25),(-5,.5),(24-1,.5),(36+2,.5)]
padrao4_1=[(0,.25),(-5,.25),(-9,.5),(-12,.5),(-5,.5)]
pp_+=padrao1+padrao2+padrao3+padrao4_1

#Ciclo3
padrao1=[(i,.5) for i in [0,12,24,36]]
padrao2=[(0,.25),(0,.25),(0,.5),(24+2,.5),(36+9,.25),(36+12,.25)]
padrao3=[(11,.25),(-1,.25),(7,.5),(24+7,.5),(36+2,.25),(36-1,.25)]
padrao5=[(0,.25),(-3,.25),(4,.5),(9,.5),(-5,.5)]
padrao5b=[(0,.25),(-3,.25),(4,.25),(9,.25),(-5,.25),(0,.25),(-3,.25),(4,.25)]
padrao4=[(0,.25),(2,.25),(2,.5),(24+2,.5),(36+9,.25),(-36+9,.25)]
padrao3_=[(11,.25),(-1,.25),(7,.25),(24+7,.5),(36+2,.5),(36-1,.25)]
padrao=   [(i,.5) for i in [0,0,-1,0]]
padrao1_= [(i,.5) for i in [0,4,7,12]]
padrao2_= [(i,.5) for i in [0,-5,-9,-12]]
padrao3__=[(i,.5) for i in [-9,-5,-3,-1]]
#padraoF=  [(i,.5) for i in [-12]]
padraoF=  [(-12,.75),(-12,.25),(-12,.75),(-12-5,.25),(-12-12,2)]
pp_+=padrao1+padrao2+padrao3+padrao5+padrao5b+padrao4+padrao3_+padrao+padrao1_+padrao2_+padrao3__+padraoF

som_=n.hstack([adsr(v(f=f0*2**(pi[0]/12.),d=pi[1],tab=Q_i)) for pi in pp_])

som=n.hstack((som,som_))


# ve a diferenca enre s_ e som, adiciona zeros em s_ e soma ambos.
# coloca s__ no comeco de tudo. Grava wav
s=n.hstack((s__,   n.hstack((  s_,n.zeros(len(som)-len(s_))   ))+som ))
###s= n.hstack((  s_,n.zeros(len(som)-len(s_))   ))+som 

s=((s-s.min())/(s.max()-s.min()))*2.-1.
orig=n.copy(s)
s = n.int16(s * float(2**15-1))
w.write("lunhaniAlpha.wav",f_a, s) # escrita do som

import os, random
palavras=["lunhani","javascript","vivace","coffescript","cravelho","cravelhoooooo","craveeeeeeelho","craaaaaaaaavelho","craaaaaaaaavvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvkkkkkkkkkkklllllllllllllllhhhhhhhhhhhhhhtttttttttttpbsadoijasdfoijsaf iajsdfoasijdfasdiofjolvelhovelhooooooooooooooooooadsiajisdasdyyyyyyo yo yo yooooooooooooo"]
vozes="f3,f2,f1,f5,m5,m1,m3".split(",")
for palavra in palavras:
    os.system("espeak -vpt-pt+%s -w%s.wav '%s'"%(random.sample(vozes,1)[0],palavra[:20],palavra))
ff=w.read("javascript.wav")[1]
ff_=n.fft.fft(ff)
s=ff2=n.fft.ifft( n.hstack((ff_,n.zeros(len(ff_)) ))  ).real
js=((s-s.min())/(s.max()-s.min()))*2.-1.

ff=w.read("lunhaniFoo.wav")[1]
ff_=n.fft.fft(ff)
s=ff2=n.fft.ifft( n.hstack((ff_,n.zeros(len(ff_)) ))  ).real
lh=((s-s.min())/(s.max()-s.min()))*2.-1.

ff=w.read("vivace.wav")[1]
ff_=n.fft.fft(ff)
s=ff2=n.fft.ifft( n.hstack((ff_,n.zeros(len(ff_)) ))  ).real
viv=((s-s.min())/(s.max()-s.min()))*2.-1.

ff=w.read("coffeescript.wav")[1]
ff_=n.fft.fft(ff)
s=ff2=n.fft.ifft( n.hstack((ff_,n.zeros(len(ff_)) ))  ).real
coff=((s-s.min())/(s.max()-s.min()))*2.-1.

ff=w.read("cravelho.wav")[1]
ff_=n.fft.fft(ff)
s=ff2=n.fft.ifft( n.hstack((ff_,n.zeros(len(ff_)) ))  ).real
crav=((s-s.min())/(s.max()-s.min()))*2.-1.

ff=w.read("craaaaaaaaavvvvvvvvv.wav")[1]; ff_=n.fft.fft(ff)
s=ff2=n.fft.ifft( n.hstack((ff_,n.zeros(len(ff_)) ))  ).real
crav_=((s-s.min())/(s.max()-s.min()))*2.-1.

ff=w.read("cravelhoooooo.wav")[1]; ff_=n.fft.fft(ff)
s=ff2=n.fft.ifft( n.hstack((ff_,n.zeros(len(ff_)) ))  ).real
cravo=((s-s.min())/(s.max()-s.min()))*2.-1.

ff=w.read("craveeeeeeelho.wav")[1]; ff_=n.fft.fft(ff)
s=ff2=n.fft.ifft( n.hstack((ff_,n.zeros(len(ff_)) ))  ).real
crave=((s-s.min())/(s.max()-s.min()))*2.-1.

ff=w.read("craaaaaaaaavvvvvvvvv.wav")[1]; ff_=n.fft.fft(ff)
s=ff2=n.fft.ifft( n.hstack((ff_,n.zeros(len(ff_)) ))  ).real
crava=((s-s.min())/(s.max()-s.min()))*2.-1.

orig[10:10+len(viv)]+=viv
orig[110:110+len(coff)]+=coff
orig[1100:1100+len(lh)]+=lh
orig[11000:11000+len(js)]+=js
T=44110*0.5

TT=T*4;TK=coff
orig[TT:TT+len(TK)]+=TK
TT=T*12;TK=lh
orig[TT:TT+len(TK)]+=TK
TT=T*16;TK=viv
orig[TT:TT+len(TK)]+=TK
TT=T*20;TK=js
orig[TT:TT+len(TK)]+=TK
TT=T*21;TK=coff
orig[TT:TT+len(TK)]+=TK
TT=T*23;TK=lh
orig[TT:TT+len(TK)]+=TK
TT=T*24;TK=viv
orig[TT:TT+len(TK)]+=TK
TT=T*24.75;TK=js
orig[TT:TT+len(TK)]+=TK
TT=T*25;TK=coff
orig[TT:TT+len(TK)]+=TK


TT=T*28;TK=crav
orig[TT:TT+len(TK)]+=TK
TT=T*28.5;TK=crav
orig[TT:TT+len(TK)]+=TK
TT=T*29;TK=crav[::-1]
orig[TT:TT+len(TK)]+=TK
TT=T*30;TK=crav[::-1]
orig[TT:TT+len(TK)]+=TK
TT=T*31;TK=crav
orig[TT:TT+len(TK)]+=TK


TT=T*36;TK=crave
orig[TT:TT+len(TK)]+=TK
TT=T*36.5;TK=crave
orig[TT:TT+len(TK)]+=TK
TT=T*37;TK=crav[::-1][:len(crav)/2]
orig[TT:TT+len(TK)]+=TK
TT=T*38;TK=crav[::-1][len(crav)/2:]
orig[TT:TT+len(TK)]+=TK
TT=T*38.5;TK=crava[::-1]
orig[TT:TT+len(TK)]+=TK
TT=T*39;TK=cravo+cravo[::-1]
orig[TT:TT+len(TK)]+=TK


TT=T*40;TK=crav[len(crav)/2:]
orig[TT:TT+len(TK)]+=TK
TT=T*41;TK=crav[:len(crav)/2]
orig[TT:TT+len(TK)]+=TK
TT=T*42;TK=crav[::-1][:len(crav)/2]
orig[TT:TT+len(TK)]+=TK
TT=T*43;TK=crav[::-1][len(crav)/2:]
orig[TT:TT+len(TK)]+=TK

TT=T*40;TK=crave
orig[TT:TT+len(TK)]+=TK
TT=T*41;TK=crava
orig[TT:TT+len(TK)]+=TK
TT=T*42;TK=cravo[::-1]
orig[TT:TT+len(TK)]+=TK
TT=T*43;TK=crav_[::-1]
orig[TT:TT+len(TK)]+=TK



TT=T*48;TK=crav[len(crav)/2:]
orig[TT:TT+len(TK)]+=TK
TT=T*49;TK=crav[:len(crav)/2]
orig[TT:TT+len(TK)]+=TK
TT=T*50;TK=crav[::-1][:len(crav)/2]
orig[TT:TT+len(TK)]+=TK
TT=T*51;TK=crav[::-1][len(crav)/2:]
orig[TT:TT+len(TK)]+=TK

TT=T*52;TK=crav_
orig[TT:TT+len(TK)]+=TK
TT=T*56;TK=crav_[::-1]
orig[TT:TT+len(TK)]+=TK


TT=T*4;TK=coff[::-1]
TT+=T*80
orig[TT:TT+len(TK)]+=TK
TT=T*12;TK=lh[::-1]
TT+=T*80
orig[TT:TT+len(TK)]+=TK
TT=T*16;TK=viv[::-1]
TT+=T*80
orig[TT:TT+len(TK)]+=TK
TT=T*20;TK=js[::-1]
TT+=T*80
orig[TT:TT+len(TK)]+=TK
TT=T*21;TK=coff
TT+=T*80
orig[TT:TT+len(TK)]+=TK
TT=T*23;TK=lh[::-1]
TT+=T*80
orig[TT:TT+len(TK)]+=TK
TT=T*24;TK=viv
TT+=T*80
orig[TT:TT+len(TK)]+=TK
TT=T*24.75;TK=js[::-1]
TT+=T*80
orig[TT:TT+len(TK)]+=TK
TT=T*25;TK=coff[::-1]
TT+=T*80
orig[TT:TT+len(TK)]+=TK
TT=T*25.3;TK=coff
TT+=T*80
orig[TT:TT+len(TK)]+=TK
s=orig
s=((s-s.min())/(s.max()-s.min()))*2.-1.
s = n.int16(s * float(2**15-1))
w.write("lunhani.wav",f_a, s) # escrita do som
