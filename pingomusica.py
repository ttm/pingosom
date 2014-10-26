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
fooXY=n.linspace(0,2*n.pi,Lt,endpoint=False)
S_i=n.sin(fooXY) # um período da senóide com T amostras

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

    Gammav_i=n.floor((ii*fv*Lv)/f_a) # índices para a LUT
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

def N(arr,xx=1.):
    r=arr
    r = (((r-r.min())/(r.max()-r.min()))*2-1)*xx
    return n.int16(r * float(2**15-1))
def NN(arr):
    return 2*((arr-arr.min())/(arr.max()-arr.min()))-1
vozes="f3,f2,f1,f5,m5,m1,m3".split(",")
def fala(frase="Semicondutor livre",ss=160):
    arq=frase.split()[0]
    #os.system("espeak -vpt-pt+%s -w%s.wav -g110 -p99 -s110 -b=1 '%s'"%(random.sample(vozes,1)[0],arq,frase))
    os.system(u"espeak -vpt-pt+%s -w%s.wav -p99 -b=1 '%s' -s%i"%(random.sample(vozes,1)[0],arq,frase,ss))
    #os.system(u"espeak "+ frase +(u" -vpt-pt+%s -w%s.wav -p99 -b=1 -s%i"%(random.sample(vozes,1)[0],arq,ss)))
    #os.system("espeak -vpt-pt+%s -w%s.wav -g110 -p99 -s130 -b=1 '%s'"%(random.sample(vozes,1)[0],arq,frase))
    ff=w.read("%s.wav"%(arq,))[1]
    ff_=n.fft.fft(ff)
    s=ff2=n.fft.ifft( n.hstack((ff_,n.zeros(len(ff_)) ))  ).real
    sc_aud=((s-s.min())/(s.max()-s.min()))*2.-1.
    return sc_aud*10

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
rb=r
r = n.int16(r * float(2**15-1))
w.write('branco.wav', f_a, r)

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
rm=r
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
r = ((r-r.min())/(r.max()-r.min()))*2-1
ra=r
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
r = ((r-r.min())/(r.max()-r.min()))*2-1
rv=r
r = n.int16(r * float(2**15-1))
w.write('violeta.wav', f_a, r)


### 2.51 Ruído rosa
# a cada oitava, perde-se 3dB
fator = 10.**(-3/20.)
alphai = fator**(n.log2(fi[i0:]/f0))

c = n.copy(coefs)
c[i0:] = coefs[i0:]*alphai
# real par, imaginaria impar
c[Lambda/2+1:] = n.real(c[1:Lambda/2])[::-1] - 1j * \
    n.imag(c[1:Lambda/2])[::-1]

ruido = n.fft.ifft(c)
r = n.real(ruido)
r = ((r-r.min())/(r.max()-r.min()))*2-1
rr=r
r = n.int16(r * float(2**15-1))
w.write('rosa.wav', f_a, r)


fator = 10.**(-9/20.)
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
rp=r
r = n.int16(r * float(2**15-1))
w.write('preto.wav', f_a, r)



w.write('respira92.wav', f_a, N(H((
   adsr(ra[   :int(f_a*2.5)],S=-.5,A=1160.),
   adsr(rv[:int(f_a*2.5)],S=-.5,A=1160.),
   adsr(ra[   :int(f_a*2.5)],S=-.5,A=1160.),
   adsr(rv[:int(f_a*2.5)],S=-.5,A=1160.),
   adsr(ra[   :int(f_a*2.5)],S=-.5,A=1160.),
   adsr(rv[:int(f_a*2.5)],S=-.5,A=1160.),
                                ))))


t1=H((
           adsr(v(f=440.*(2.**(0/12.))  ,fv=0., nu=0.,  d=.5, tab=Tr_i), S=-9.,R=10.,A=5., D=40.),
           adsr(v(f=440.*(2.**(4./12.)),fv=0., nu=0.,  d=.5,  tab=Tr_i), S=-9.,R=5. ,A=5., D=40.),
           adsr(v(f=440.*(2.**(8/12.))  ,fv=0., nu=0.,  d=.5, tab=Tr_i), S=-9.,R=10.,A=5., D=40.),

           adsr(v(f=440.*(2.**(4/12.))  ,fv=0., nu=0.,  d=.5, tab=Tr_i), S=-9.,R=10.,A=5., D=40.),
           adsr(v(f=440.*(2.**(0/12.)) ,fv=0., nu=0.,  d=.5,  tab=Tr_i), S=-9.,R=10.,A=5., D=40.),
           adsr(v(f=440.*(2.**(8./12.)),fv=0., nu=0.,  d=.5,  tab=Tr_i), S=-9.,R=5. ,A=5., D=40.),

           adsr(v(f=440.*(2.**(4/12.))  ,fv=0., nu=0.,  d=.5, tab=Tr_i), S=-9.,R=10.,A=5., D=40.),
           adsr(v(f=440.*(2.**(8/12.)) ,fv=0., nu=0.,  d=.5,  tab=Tr_i), S=-9.,R=10.,A=5., D=40.),
           adsr(v(f=440.*(2.**(0./12.)),fv=0., nu=0.,  d=.5,  tab=Tr_i), S=-9.,R=5. ,A=5., D=40.),

           adsr(v(f=440.*(2.**(8/12.))  ,fv=0., nu=0.,  d=.5, tab=Tr_i), S=-9.,R=10.,A=5., D=40.),
           adsr(v(f=440.*(2.**(4/12.)) ,fv=0., nu=0.,  d=.5,  tab=Tr_i), S=-9.,R=10.,A=5., D=40.),
           adsr(v(f=440.*(2.**(0./12.)),fv=0., nu=0.,  d=.5,  tab=Tr_i), S=-9.,R=5. ,A=5., D=40.),

           adsr(v(f=440.*(2.**(8/12.))  ,fv=0., nu=0.,  d=.5, tab=Tr_i), S=-9.,R=10.,A=5., D=40.),
           adsr(v(f=440.*(2.**(0/12.)) ,fv=0., nu=0.,  d=.5,  tab=Tr_i), S=-9.,R=10.,A=5., D=40.),
           adsr(v(f=440.*(2.**(4./12.)),fv=0., nu=0.,  d=.5,  tab=Tr_i), S=-9.,R=5. ,A=5., D=40.),

           adsr(v(f=440.*(2.**(0/12.))  ,fv=0., nu=0.,  d=.5, tab=Tr_i), S=-9.,R=10.,A=5., D=40.),
           adsr(v(f=440.*(2.**(8/12.)) ,fv=0., nu=0.,  d=.5,  tab=Tr_i), S=-9.,R=10.,A=5., D=40.),
           adsr(v(f=440.*(2.**(4./12.)),fv=0., nu=0.,  d=.5,  tab=Tr_i), S=-9.,R=5. ,A=5., D=40.),
))

t2=H((
           adsr(v(f=110.*(2.**(0/12.))  ,fv=8., nu=0.2, d=1.5, tab=Q_i), S=-9.,R=10.,A=5., D=40.),
           adsr(v(f=110.*(2.**(4./12.)), fv=8., nu=0.2,  d=1.5, tab=Q_i), S=-9.,R=5. ,A=5., D=40.),
           adsr(v(f=110.*(2.**(8/12.))  ,fv=8., nu=0.2, d=1.5, tab=Q_i), S=-9.,R=10.,A=5., D=40.),

           adsr(v(f=110.*(2.**(4/12.))  ,fv=8., nu=0.2, d=1.5, tab=Q_i), S=-9.,R=10.,A=5., D=40.),
           adsr(v(f=110.*(2.**(0./12.)), fv=8., nu=0.2,  d=1.5, tab=Q_i), S=-9.,R=5. ,A=5., D=40.),
           adsr(v(f=110.*(2.**(8/12.))  ,fv=8., nu=0.2, d=1.5, tab=Q_i), S=-9.,R=10.,A=5., D=40.),

           ## novo ciclo
           adsr(v(f=110.*(2.**(4/12.))  ,fv=8., nu=0.2, d=1.5, tab=Q_i), S=-9.,R=10.,A=5., D=40.),
           adsr(v(f=110.*(2.**(8./12.)), fv=8., nu=0.2,  d=1.5, tab=Q_i), S=-9.,R=5. ,A=5., D=40.),
           adsr(v(f=110.*(2.**(0/12.))  ,fv=8., nu=0.2, d=1.5, tab=Q_i), S=-9.,R=10.,A=5., D=40.),

           adsr(v(f=110.*(2.**(8/12.))  ,fv=8., nu=0.2, d=1.5, tab=Q_i), S=-9.,R=10.,A=5., D=40.),
           adsr(v(f=110.*(2.**(4./12.)), fv=8., nu=0.2,  d=1.5, tab=Q_i), S=-9.,R=5. ,A=5., D=40.),
           adsr(v(f=110.*(2.**(0/12.))  ,fv=8., nu=0.2, d=1.5, tab=Q_i), S=-9.,R=10.,A=5., D=40.),

           ## novo ciclo
           adsr(v(f=110.*(2.**(8/12.))  ,fv=8., nu=0.2, d=1.5, tab=Q_i), S=-9.,R=10.,A=5., D=40.),
           adsr(v(f=110.*(2.**(0./12.)), fv=8., nu=0.2,  d=1.5, tab=Q_i), S=-9.,R=5. ,A=5., D=40.),
           adsr(v(f=110.*(2.**(4/12.))  ,fv=8., nu=0.2, d=1.5, tab=Q_i), S=-9.,R=10.,A=5., D=40.),

           adsr(v(f=110.*(2.**(0/12.))  ,fv=8., nu=0.2, d=1.5, tab=Q_i), S=-9.,R=10.,A=5., D=40.),
           adsr(v(f=110.*(2.**(8./12.)), fv=8., nu=0.2,  d=1.5, tab=Q_i), S=-9.,R=5. ,A=5., D=40.),
           adsr(v(f=110.*(2.**(4/12.))  ,fv=8., nu=0.2, d=1.5, tab=Q_i), S=-9.,R=10.,A=5., D=40.),
))


w.write('musica1.wav', f_a, N(H((t1,t1,t1))+t2))
