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


#w.write('respira.wav', f_a, N(H((
#                                 rr[:int(f_a*.5)],
#                                 rm[:int(f_a*.5)],
#                                 rr[:int(f_a*.5)],
#                                 rm[:int(f_a*.5)],
#                                 rr[:int(f_a*.5)],
#                                 rm[:int(f_a*.5)],
#                                ))))
#
#w.write('respira2.wav', f_a, N(H((
#                                 rp[:int(f_a*.5)],
#                                 rm[:int(f_a*.5)],
#                                 rp[:int(f_a*.5)],
#                                 rm[:int(f_a*.5)],
#                                 rp[:int(f_a*.5)],
#                                 rm[:int(f_a*.5)],
#                                ))))
#
#
#w.write('respira3.wav', f_a, N(H((
#                                 adsr(rr[:int(f_a*.5)],S=-.5,A=360.),
#                                 5.*adsr(rm[:int(f_a*.5)],S=-.5,A=360.),
#                                 adsr(rr[:int(f_a*.5)],S=-.5,A=360.),
#                                 5.*adsr(rm[:int(f_a*.5)],S=-.5,A=360.),
#                                 adsr(rr[:int(f_a*.5)],S=-.5,A=360.),
#                                 5.*adsr(rm[:int(f_a*.5)],S=-.5,A=360.),
#                                ))))
#
#
#w.write('respira4.wav', f_a, N(H((
#                                 adsr(rr[:int(f_a*.5)],S=-.5,A=360.),
#                                 adsr(rb[:int(f_a*.5)],S=-.5,A=360.),
#                                 adsr(rr[:int(f_a*.5)],S=-.5,A=360.),
#                                 adsr(rb[:int(f_a*.5)],S=-.5,A=360.),
#                                 adsr(rr[:int(f_a*.5)],S=-.5,A=360.),
#                                 adsr(rb[:int(f_a*.5)],S=-.5,A=360.),
#                                ))))
#
#
#w.write('respira5.wav', f_a, N(H((
#                                 adsr(rr[:int(f_a*.5)],S=-.5,A=360.),
#                                 adsr(rv[:int(f_a*.5)],S=-.5,A=360.),
#                                 adsr(rr[:int(f_a*.5)],S=-.5,A=360.),
#                                 adsr(rv[:int(f_a*.5)],S=-.5,A=360.),
#                                 adsr(rr[:int(f_a*.5)],S=-.5,A=360.),
#                                 adsr(rv[:int(f_a*.5)],S=-.5,A=360.),
#                                ))))
#
#
#w.write('respira6.wav', f_a, N(H((
#                                 adsr(rr[:int(f_a*.2)],S=-.5,A=160.,R=10.),
#                                 adsr(rv[:int(f_a*.2)],S=-.5,A=160.,R=10.),
#                                 adsr(rr[:int(f_a*.2)],S=-.5,A=160.,R=10.),
#                                 adsr(rv[:int(f_a*.2)],S=-.5,A=160.,R=10.),
#                                 adsr(rr[:int(f_a*.2)],S=-.5,A=160.,R=10.),
#                                 adsr(rv[:int(f_a*.2)],S=-.5,A=160.,R=10.),
#                                 adsr(rr[:int(f_a*.2)],S=-.5,A=160.,R=10.),
#                                 adsr(rv[:int(f_a*.2)],S=-.5,A=160.,R=10.),
#                                 adsr(rr[:int(f_a*.2)],S=-.5,A=160.,R=10.),
#                                 adsr(rv[:int(f_a*.2)],S=-.5,A=160.,R=10.),
#                                 adsr(rr[:int(f_a*.2)],S=-.5,A=160.,R=10.),
#                                 adsr(rv[:int(f_a*.2)],S=-.5,A=160.,R=10.),
#                                 adsr(rr[:int(f_a*.2)],S=-.5,A=160.,R=10.),
#                                 adsr(rv[:int(f_a*.2)],S=-.5,A=160.,R=10.),
#                                 adsr(rr[:int(f_a*.2)],S=-.5,A=160.,R=10.),
#                                 adsr(rv[:int(f_a*.2)],S=-.5,A=160.,R=10.),
#                                ))))
#
#
#f0=110.
#s=n.zeros(4*f_a)
#kk=(2*n.pi/10)*2. # uma volta
#aa=20. # 10. dB
#for i in xrange(10): # 10 harmonicas
#    ff=f0*(1+i)
#    n_oitavas=n.log2(ff/f0)
#    a_=10.**((n_oitavas*(-25.+aa*n.cos(kk*i))/20.))
#    s+=v(f=ff,d=4.,nu=0.)*a_
#
#w.write('pisca.wav',  f_a, N(H((
#                              s[:f_a/2], n.zeros(f_a/2),
#                              s[:f_a/2], n.zeros(f_a/2),
#                              s[:f_a/2], n.zeros(f_a/2),
#                              s[:f_a/2], n.zeros(f_a/2),
#                           ))))
#
#
#
#f0=1100.
#s=n.zeros(4*f_a)
#kk=(2*n.pi/10)*2. # uma volta
#aa=20. # 10. dB
#for i in xrange(10): # 10 harmonicas
#    ff=f0*(1+i)
#    n_oitavas=n.log2(ff/f0)
#    a_=10.**((n_oitavas*(-25.+aa*n.cos(kk*i))/20.))
#    s+=v(f=ff,d=4.,nu=0.)*a_
#
#w.write('pisca2.wav',  f_a, N(H((
#                              s[:f_a/2], n.zeros(f_a/2),
#                              s[:f_a/2], n.zeros(f_a/2),
#                              s[:f_a/2], n.zeros(f_a/2),
#                              s[:f_a/2], n.zeros(f_a/2),
#                           ))))
#
#
#
#f0=11000.
#s=n.zeros(4*f_a)
#kk=(2*n.pi/10)*2. # uma volta
#aa=20. # 10. dB
#for i in xrange(10): # 10 harmonicas
#    ff=f0*(1+i)
#    n_oitavas=n.log2(ff/f0)
#    a_=10.**((n_oitavas*(-25.+aa*n.cos(kk*i))/20.))
#    s+=v(f=ff,d=4.,nu=0.)*a_
#
#w.write('pisca3.wav',  f_a, N(H((
#                              s[:f_a/2], n.zeros(f_a/2),
#                              s[:f_a/2], n.zeros(f_a/2),
#                              s[:f_a/2], n.zeros(f_a/2),
#                              s[:f_a/2], n.zeros(f_a/2),
#                           ))))
#
#
#
#f0=410.
#s=n.zeros(4*f_a)
#kk=(2*n.pi/10)*2. # uma volta
#aa=20. # 10. dB
#for i in xrange(10): # 10 harmonicas
#    ff=f0*(1+i)
#    n_oitavas=n.log2(ff/f0)
#    a_=10.**((n_oitavas*(-25.+aa*n.cos(kk*i))/20.))
#    s+=adsr(v(f=ff,d=4.,nu=0.)*a_,S=-5.)
#
#w.write('pisca4.wav',  f_a, N(H((
#                              s[:f_a/2], n.zeros(f_a/2),
#                              s[:f_a/2], n.zeros(f_a/2),
#                              s[:f_a/2], n.zeros(f_a/2),
#                              s[:f_a/2], n.zeros(f_a/2),
#                           ))))

##### PISCA TTMPPC
#f0=110.
#s=n.zeros(4*f_a)
#kk=(2*n.pi/10)*2. # uma volta
#aa=20. # 10. dB
#for i in xrange(10): # 10 harmonicas
#    ff=f0*(1+i)
#    n_oitavas=n.log2(ff/f0)
#    a_=10.**((n_oitavas*(-25.+aa*n.cos(kk*i))/20.))
#    s+=v(f=ff,d=4.,nu=0.)*a_
#
#w.write('pisca_.wav',  f_a, N(H((
#                              s[:f_a/8], n.zeros(f_a/2),
#                           ))))
#
#
#
#f0=1100.
#s=n.zeros(4*f_a)
#kk=(2*n.pi/10)*2. # uma volta
#aa=20. # 10. dB
#for i in xrange(10): # 10 harmonicas
#    ff=f0*(1+i)
#    n_oitavas=n.log2(ff/f0)
#    a_=10.**((n_oitavas*(-25.+aa*n.cos(kk*i))/20.))
#    s+=v(f=ff,d=4.,nu=0.)*a_
#
#w.write('pisca2_.wav',  f_a, N(H((
#                              s[:f_a/8], n.zeros(f_a/2),
#                           ))))
#
#
#
#f0=11000.
#s=n.zeros(4*f_a)
#kk=(2*n.pi/10)*2. # uma volta
#aa=20. # 10. dB
#for i in xrange(10): # 10 harmonicas
#    ff=f0*(1+i)
#    n_oitavas=n.log2(ff/f0)
#    a_=10.**((n_oitavas*(-25.+aa*n.cos(kk*i))/20.))
#    s+=v(f=ff,d=4.,nu=0.)*a_
#
#w.write('pisca3_.wav',  f_a, N(H((
#                              s[:f_a/8], n.zeros(f_a/2),
#                           ))))
#
#
#
#f0=410.
#s=n.zeros(4*f_a)
#kk=(2*n.pi/10)*2. # uma volta
#aa=20. # 10. dB
#for i in xrange(10): # 10 harmonicas
#    ff=f0*(1+i)
#    n_oitavas=n.log2(ff/f0)
#    a_=10.**((n_oitavas*(-25.+aa*n.cos(kk*i))/20.))
#    s+=adsr(v(f=ff,d=4.,nu=0.)*a_,S=-5.)
#
#w.write('pisca4_.wav',  f_a, N(H((
#                              s[:f_a/8], n.zeros(f_a/2),
#                           ))))
#

##### END TTMPPC

w.write('comendo6.wav',  f_a, N(fala("O melhor que voce faz com a sua boca, eh servir de toca, para outra cabessa. Nao que voce meressa, esta oportunidade, que vem com a idade, de se curtir em mim.",ss=3500)))

w.write('comendo7.wav',  f_a, N(fala("Diga aonde voce vai, que eu vou varrendo, diga aonda voce vai, que eu vou varrendo. Vou varrendo, vou varrendo vou varrendo. Vou varrendo, vou varrendo, vou varrendo.",ss=3500)))


#
#
#w.write('comendo.wav',  f_a, N(fala("mahnamnahamhahamnahamhanhamnanhnahamha")))
#w.write('comendo2.wav',  f_a, N(fala("manamnaamaamnaamanamnannaama")))
#w.write('comendo3.wav',  f_a, N(fala("mnmnmmnmnmnnnm")))
#w.write('comendo4.wav',  f_a, N(fala("mnmnmm nmnm nn nmnmnmn")))
#w.write('comendo5.wav',  f_a, N(fala("mnhmnhmm nhmhnm nn nhmhnmhnhmn")))
#
#
#w.write('chorando_.wav',  f_a, N(fala("bbbbuaaa bbbbbuaaa bbbbuaaa bbbuaaa")))
#
#
#w.write('chorando_2.wav',  f_a, N(fala("buaaa bbuaaa buaaa buaaa")))
#
#
#
#w.write('chorando_3.wav',  f_a, N(fala("buaaa nheee ee ee nheeee e eeeee bbuaaa buaaa nheeeee eee eeeee buaaa")))
#
#
#w.write('chorando_4.wav',  f_a, N(fala("buaaa nheee ee hhh hhh hhh ee nheeehhhh h hh hhe e eeeee bbuhhh h hh haaa buaaa nhhhh hhh eeeee eee hhhhhh h heeeee buaaa")))
#

w.write('coma.wav',  f_a, N(H((
       v(f=1000.*(3./2),nu=0.,d=0.5), n.zeros(int(f_a*0.5)),
       v(f=1000.*(3./2),nu=0.,d=0.5), n.zeros(int(f_a*0.5)),
       v(f=1000.*(3./2),nu=0.,d=0.5), n.zeros(int(f_a*0.5)),
       v(f=1000.*(3./2),nu=0.,d=0.5), n.zeros(int(f_a*0.5)),
                           )),.3))

w.write('coma2.wav',  f_a, N(H((
       v(f=1000.*2.*(3./2),nu=0.,d=0.5), n.zeros(int(f_a*0.5)),
       v(f=1000.*2.*(3./2),nu=0.,d=0.5), n.zeros(int(f_a*0.5)),
       v(f=1000.*2.*(3./2),nu=0.,d=0.5), n.zeros(int(f_a*0.5)),
       v(f=1000.*2.*(3./2),nu=0.,d=0.5), n.zeros(int(f_a*0.5)),
                           )),.3))

w.write('coma3.wav',  f_a, N(H((
       v(f=1000.*2.*(3./2),nu=0.,d=0.1), v(f=1000.*2.*(3./2),nu=0.,d=0.3), n.zeros(int(f_a*1.5)),
       v(f=1000.*2.*(3./2),nu=0.,d=0.1), v(f=1000.*2.*(3./2),nu=0.,d=0.3), n.zeros(int(f_a*1.5)),
       v(f=1000.*2.*(3./2),nu=0.,d=0.1), v(f=1000.*2.*(3./2),nu=0.,d=0.3), n.zeros(int(f_a*1.5)),
       v(f=1000.*2.*(3./2),nu=0.,d=0.1), v(f=1000.*2.*(3./2),nu=0.,d=0.3), n.zeros(int(f_a*1.5)),
                           )),.3))


w.write('coma4.wav',  f_a, N(H((
       v(f=1000.*2.*(3./2),nu=0.,d=0.1), n.zeros(int(f_a*0.1)), v(f=1000.*2.*(3./2),nu=0.,d=0.3), n.zeros(int(f_a*1.5)),
       v(f=1000.*2.*(3./2),nu=0.,d=0.1), n.zeros(int(f_a*0.1)), v(f=1000.*2.*(3./2),nu=0.,d=0.3), n.zeros(int(f_a*1.5)),
       v(f=1000.*2.*(3./2),nu=0.,d=0.1), n.zeros(int(f_a*0.1)), v(f=1000.*2.*(3./2),nu=0.,d=0.3), n.zeros(int(f_a*1.5)),
       v(f=1000.*2.*(3./2),nu=0.,d=0.1), n.zeros(int(f_a*0.1)), v(f=1000.*2.*(3./2),nu=0.,d=0.3), n.zeros(int(f_a*1.5)),
                           )),.3))


w.write('coma5.wav',  f_a, N(H((
       v(f=1000.*2.*(3./2),nu=0.,d=0.1), n.zeros(int(f_a*1.5)),
       v(f=1000.*2.*(3./2),nu=0.,d=0.1), n.zeros(int(f_a*1.5)),
       v(f=1000.*2.*(3./2),nu=0.,d=0.1), n.zeros(int(f_a*1.5)),
       v(f=1000.*2.*(3./2),nu=0.,d=0.1), n.zeros(int(f_a*1.5)),
                           )),.3))

w.write('coma6.wav',  f_a, N(H((
       v(f=1000.*2.*(3./2),nu=0.,d=0.1), n.zeros(int(f_a*2.5)),
       v(f=1000.*2.*(3./2),nu=0.,d=0.1), n.zeros(int(f_a*2.5)),
       v(f=1000.*2.*(3./2),nu=0.,d=0.1), n.zeros(int(f_a*2.5)),
       v(f=1000.*2.*(3./2),nu=0.,d=0.1), n.zeros(int(f_a*2.5)),
                           )),.3))

w.write('coma7.wav',  f_a, N(H((
       v(f=1000.*2.*(3./2),nu=0.,d=0.1), n.zeros(int(f_a*3.5)),
       v(f=1000.*2.*(3./2),nu=0.,d=0.1), n.zeros(int(f_a*3.5)),
       v(f=1000.*2.*(3./2),nu=0.,d=0.1), n.zeros(int(f_a*3.5)),
       v(f=1000.*2.*(3./2),nu=0.,d=0.1), n.zeros(int(f_a*3.5)),
                           )),.3))

w.write('coma8.wav',  f_a, N(H((
       v(f=1000.*2.*(3./2),nu=2.,d=0.1,tab=Tr_i), n.zeros(int(f_a*3.5)),
       v(f=1000.*2.*(3./2),nu=2.,d=0.1,tab=Tr_i), n.zeros(int(f_a*3.5)),
       v(f=1000.*2.*(3./2),nu=2.,d=0.1,tab=Tr_i), n.zeros(int(f_a*3.5)),
       v(f=1000.*2.*(3./2),nu=2.,d=0.1,tab=Tr_i), n.zeros(int(f_a*3.5)),
                           )),.3))



w.write('respira7.wav', f_a, N(H((
   adsr(rr[   :int(f_a*1.5)],S=-.5,A=360.),
   5.*adsr(rm[:int(f_a*1.5)],S=-.5,A=360.),
   adsr(rr[   :int(f_a*1.5)],S=-.5,A=360.),
   5.*adsr(rm[:int(f_a*1.5)],S=-.5,A=360.),
   adsr(rr[   :int(f_a*1.5)],S=-.5,A=360.),
   5.*adsr(rm[:int(f_a*1.5)],S=-.5,A=360.),
                                ))))


w.write('respira8.wav', f_a, N(H((
   adsr(rr[   :int(f_a*2.5)],S=-.5,A=360.),
   5.*adsr(rm[:int(f_a*2.5)],S=-.5,A=360.),
   adsr(rr[   :int(f_a*2.5)],S=-.5,A=360.),
   5.*adsr(rm[:int(f_a*2.5)],S=-.5,A=360.),
   adsr(rr[   :int(f_a*2.5)],S=-.5,A=360.),
   5.*adsr(rm[:int(f_a*2.5)],S=-.5,A=360.),
                                ))))

w.write('respira9.wav', f_a, N(H((
   adsr(rr[   :int(f_a*2.5)],S=-.5,A=1160.),
   adsr(rb[:int(f_a*2.5)],S=-.5,A=1160.),
   adsr(rr[   :int(f_a*2.5)],S=-.5,A=1160.),
   adsr(rb[:int(f_a*2.5)],S=-.5,A=1160.),
   adsr(rr[   :int(f_a*2.5)],S=-.5,A=1160.),
   adsr(rb[:int(f_a*2.5)],S=-.5,A=1160.),
                                ))))


w.write('respira91.wav', f_a, N(H((
   adsr(ra[   :int(f_a*2.5)],S=-.5,A=1160.),
   adsr(rb[:int(f_a*2.5)],S=-.5,A=1160.),
   adsr(ra[   :int(f_a*2.5)],S=-.5,A=1160.),
   adsr(rb[:int(f_a*2.5)],S=-.5,A=1160.),
   adsr(ra[   :int(f_a*2.5)],S=-.5,A=1160.),
   adsr(rb[:int(f_a*2.5)],S=-.5,A=1160.),
                                ))))

w.write('respira92.wav', f_a, N(H((
   adsr(ra[   :int(f_a*2.5)],S=-.5,A=1160.),
   adsr(rv[:int(f_a*2.5)],S=-.5,A=1160.),
   adsr(ra[   :int(f_a*2.5)],S=-.5,A=1160.),
   adsr(rv[:int(f_a*2.5)],S=-.5,A=1160.),
   adsr(ra[   :int(f_a*2.5)],S=-.5,A=1160.),
   adsr(rv[:int(f_a*2.5)],S=-.5,A=1160.),
                                ))))


w.write('dormindo.wav', f_a, N(H((
   adsr(ra[   :int(f_a*2.5)],S=-.5,A=1160.),
   adsr(rv[:int(f_a*2.5)],S=-.5,A=1160.),
   adsr(ra[   :int(f_a*2.5)],S=-.5,A=1160.),
   adsr(rv[:int(f_a*2.5)],S=-.5,A=1160.),
   adsr(ra[   :int(f_a*2.5)],S=-.5,A=1160.),
   adsr(rv[:int(f_a*2.5)],S=-.5,A=1160.),
                                ))))

# arroto3 arroto6 arroto 9 92

w.write('dormindo2.wav', f_a, N(H((
   adsr(ra,S=-.5,A=1760.),
   adsr(rv,S=-.5,A=1760.),
   adsr(ra,S=-.5,A=1760.),
   adsr(rv,S=-.5,A=1760.),
   adsr(ra,S=-.5,A=1760.),
   adsr(rv,S=-.5,A=1760.),
                                ))))

w.write('dormindo2.wav', f_a, N(H((
   adsr(ra,S=-.5,A=1760.),
   adsr(rv,S=-.5,A=1760.),
   adsr(ra,S=-.5,A=1760.),
   adsr(rv,S=-.5,A=1760.),
   adsr(ra,S=-.5,A=1760.),
   adsr(rv,S=-.5,A=1760.),
                                ))))

ronco=H((
                         adsr(rp[:int(f_a*0.040)],A=3.,S=-3.,R=10.),
                         adsr(rp[:int(f_a*0.040)],A=3.,S=-3.,R=10.),
                         adsr(rp[:int(f_a*0.040)],A=3.,S=-3.,R=10.),
                         adsr(rp[:int(f_a*0.040)],A=3.,S=-3.,R=10.),
                         adsr(rp[:int(f_a*0.040)],A=3.,S=-3.,R=10.),
                         adsr(rp[:int(f_a*0.040)],A=3.,S=-3.,R=10.),
                         adsr(rp[:int(f_a*0.040)],A=3.,S=-3.,R=10.),
                      ))

w.write('dormindo3.wav', f_a, N(H((
                                   ronco,n.zeros(f_a),
                                   ronco,n.zeros(f_a),
                                   ronco,n.zeros(f_a),
                                   ronco,n.zeros(f_a),
                                   ronco,n.zeros(f_a),
                                 ))))

w.write('dormindo4.wav', f_a, N(H((
   adsr(ra,S=-.5,A=1760.),ronco,n.zeros(f_a),
   adsr(rv,S=-.5,A=1760.),ronco,n.zeros(f_a),
   adsr(ra,S=-.5,A=1760.),ronco,n.zeros(f_a),
   adsr(rv,S=-.5,A=1760.),ronco,n.zeros(f_a),
   adsr(ra,S=-.5,A=1760.),ronco,n.zeros(f_a),
   adsr(rv,S=-.5,A=1760.),ronco,n.zeros(f_a),
                                ))))


w.write('dormindo5.wav', f_a, N(H((
   adsr(ra,S=-.5,A=1760.),10*ronco,n.zeros(f_a),
   adsr(rv,S=-.5,A=1760.),10*ronco,n.zeros(f_a),
   adsr(ra,S=-.5,A=1760.),10*ronco,n.zeros(f_a),
   adsr(rv,S=-.5,A=1760.),10*ronco,n.zeros(f_a),
   adsr(ra,S=-.5,A=1760.),10*ronco,n.zeros(f_a),
   adsr(rv,S=-.5,A=1760.),10*ronco,n.zeros(f_a),
                                ))))


w.write('dormindo6.wav', f_a, N(H((
   adsr(ra,S=-.5,A=1760.),5*ronco,n.zeros(f_a),
   adsr(rv,S=-.5,A=1760.),5*ronco,n.zeros(f_a),
   adsr(ra,S=-.5,A=1760.),5*ronco,n.zeros(f_a),
   adsr(rv,S=-.5,A=1760.),5*ronco,n.zeros(f_a),
   adsr(ra,S=-.5,A=1760.),5*ronco,n.zeros(f_a),
   adsr(rv,S=-.5,A=1760.),5*ronco,n.zeros(f_a),
                                ))))
w.write('dormindo7.wav', f_a, N(H((
   adsr(ra,S=-.5,A=1760.)+H((n.zeros(len(ra)-len(ronco)),5*ronco)),n.zeros(f_a),
   adsr(rv,S=-.5,A=1760.)+H((n.zeros(len(ra)-len(ronco)),5*ronco)),n.zeros(f_a),
   adsr(ra,S=-.5,A=1760.)+H((n.zeros(len(ra)-len(ronco)),5*ronco)),n.zeros(f_a),
   adsr(rv,S=-.5,A=1760.)+H((n.zeros(len(ra)-len(ronco)),5*ronco)),n.zeros(f_a),
   adsr(ra,S=-.5,A=1760.)+H((n.zeros(len(ra)-len(ronco)),5*ronco)),n.zeros(f_a),
   adsr(rv,S=-.5,A=1760.)+H((n.zeros(len(ra)-len(ronco)),5*ronco)),n.zeros(f_a),
                                ))))


ronco2=H((
                         adsr(rm[:int(f_a*0.040)],A=5.,S=-3.,R=10.),
                         adsr(rm[:int(f_a*0.040)],A=5.,S=-3.,R=10.),
                         adsr(rm[:int(f_a*0.040)],A=5.,S=-3.,R=10.),
                         adsr(rm[:int(f_a*0.040)],A=5.,S=-3.,R=10.),
                         adsr(rm[:int(f_a*0.040)],A=5.,S=-3.,R=10.),
                         adsr(rm[:int(f_a*0.040)],A=5.,S=-3.,R=10.),
                         adsr(rm[:int(f_a*0.040)],A=5.,S=-3.,R=10.),
                      ))

w.write('dormindo8.wav', f_a, N(H((
   adsr(ra,S=-.5,A=1760.)+H((n.zeros(len(ra)-len(ronco)),5*ronco)),n.zeros(f_a),ronco2,
   adsr(rv,S=-.5,A=1760.)+H((n.zeros(len(ra)-len(ronco)),5*ronco)),n.zeros(f_a),ronco2,
   adsr(ra,S=-.5,A=1760.)+H((n.zeros(len(ra)-len(ronco)),5*ronco)),n.zeros(f_a),ronco2,
   adsr(rv,S=-.5,A=1760.)+H((n.zeros(len(ra)-len(ronco)),5*ronco)),n.zeros(f_a),ronco2,
   adsr(ra,S=-.5,A=1760.)+H((n.zeros(len(ra)-len(ronco)),5*ronco)),n.zeros(f_a),ronco2,
   adsr(rv,S=-.5,A=1760.)+H((n.zeros(len(ra)-len(ronco)),5*ronco)),n.zeros(f_a),ronco2,
                                ))))


w.write('dormindo9.wav', f_a, N(H((
   adsr(ra,S=-.5,A=1760.)+H((n.zeros(len(ra)-len(ronco)),5*ronco)),v(440.,nu=3.,fv=.2),ronco2,
   adsr(rv,S=-.5,A=1760.)+H((n.zeros(len(ra)-len(ronco)),5*ronco)),v(440.,nu=3.,fv=.2),ronco2,
   adsr(ra,S=-.5,A=1760.)+H((n.zeros(len(ra)-len(ronco)),5*ronco)),v(440.,nu=3.,fv=.2),ronco2,
   adsr(rv,S=-.5,A=1760.)+H((n.zeros(len(ra)-len(ronco)),5*ronco)),v(440.,nu=3.,fv=.2),ronco2,
   adsr(ra,S=-.5,A=1760.)+H((n.zeros(len(ra)-len(ronco)),5*ronco)),v(440.,nu=3.,fv=.2),ronco2,
   adsr(rv,S=-.5,A=1760.)+H((n.zeros(len(ra)-len(ronco)),5*ronco)),v(440.,nu=3.,fv=.2),ronco2,
                                ))))


w.write('dormindo91.wav', f_a, N(H((
   adsr(ra,S=-.5,A=1760.)+H((n.zeros(len(ra)-len(ronco)),5*ronco)),v(440.*2,nu=3.,fv=.2,tabv=D_i),ronco2,
   adsr(rv,S=-.5,A=1760.)+H((n.zeros(len(ra)-len(ronco)),5*ronco)),v(440.*2,nu=3.,fv=.2,tabv=D_i),ronco2,
   adsr(ra,S=-.5,A=1760.)+H((n.zeros(len(ra)-len(ronco)),5*ronco)),v(440.*2,nu=3.,fv=.2,tabv=D_i),ronco2,
   adsr(rv,S=-.5,A=1760.)+H((n.zeros(len(ra)-len(ronco)),5*ronco)),v(440.*2,nu=3.,fv=.2,tabv=D_i),ronco2,
   adsr(ra,S=-.5,A=1760.)+H((n.zeros(len(ra)-len(ronco)),5*ronco)),v(440.*2,nu=3.,fv=.2,tabv=D_i),ronco2,
   adsr(rv,S=-.5,A=1760.)+H((n.zeros(len(ra)-len(ronco)),5*ronco)),v(440.*2,nu=3.,fv=.2,tabv=D_i),ronco2,
                                ))))


w.write('dormindo92.wav', f_a, N(H((
   adsr(ra,S=-.5,A=1760.)+H((n.zeros(len(ra)-len(ronco)),5*ronco)),v(440.*2,nu=3.,fv=.2,tabv=D_i[::-1]),ronco2,
   adsr(rv,S=-.5,A=1760.)+H((n.zeros(len(ra)-len(ronco)),5*ronco)),v(440.*2,nu=3.,fv=.2,tabv=D_i[::-1]),ronco2,
   adsr(ra,S=-.5,A=1760.)+H((n.zeros(len(ra)-len(ronco)),5*ronco)),v(440.*2,nu=3.,fv=.2,tabv=D_i[::-1]),ronco2,
   adsr(rv,S=-.5,A=1760.)+H((n.zeros(len(ra)-len(ronco)),5*ronco)),v(440.*2,nu=3.,fv=.2,tabv=D_i[::-1]),ronco2,
   adsr(ra,S=-.5,A=1760.)+H((n.zeros(len(ra)-len(ronco)),5*ronco)),v(440.*2,nu=3.,fv=.2,tabv=D_i[::-1]),ronco2,
   adsr(rv,S=-.5,A=1760.)+H((n.zeros(len(ra)-len(ronco)),5*ronco)),v(440.*2,nu=3.,fv=.2,tabv=D_i[::-1]),ronco2,
                                ))))

w.write('porta_abre.wav', f_a, N(v(200,fv=1./(7*2.),d=1.0,nu=20.)))
w.write('porta_abre2.wav', f_a, N(v(800,fv=1./(7*2.),d=1.0,nu=20.)))
w.write('porta_abre3.wav', f_a, N(v(800,fv=1.,d=.5,nu=20.,tabv=D_i)))
w.write('porta_abre4.wav', f_a, N(v(1800,fv=1.,d=.5,nu=20.,tabv=D_i,tab=Tr_i)))
w.write('porta_abre5.wav', f_a, N(v(2800,fv=1.,d=.5,nu=20.,tabv=D_i,tab=Tr_i)))
w.write('porta_abre6.wav', f_a, N(v(2800,fv=1.,d=.5,nu=2.,tabv=D_i,tab=Tr_i)))
w.write('porta_abre7.wav', f_a, N(v(1800,fv=1.,d=.5,nu=20.,tabv=D_i,tab=D_i)))
w.write('porta_abre8.wav', f_a, N(v(1800,fv=1.,d=.5,nu=20.,tabv=D_i,tab=Q_i)))

w.write('porta_fecha.wav', f_a, N(v(200,fv=1./(7*2.),d=1.0,nu=20.  , tabv=S_i*-1)))
w.write('porta_fecha2.wav', f_a, N(v(800,fv=1./(7*2.),d=1.0,nu=20. , tabv=S_i*-1)))
w.write('porta_fecha3.wav', f_a, N(v(800,fv=1.,d=.5,nu=20.,tabv=D_i)))
w.write('porta_fecha4.wav', f_a, N(v(1800,fv=1.,d=.5,nu=20.,tabv=D_i,tab=Tr_i*-1)))
w.write('porta_fecha5.wav', f_a, N(v(2800,fv=1.,d=.5,nu=20.,tabv=D_i,tab=Tr_i*-1)))
w.write('porta_fecha6.wav', f_a, N(v(2800,fv=1.,d=.5,nu=2.,tabv=D_i,tab=Tr_i *-1)))
w.write('porta_fecha7.wav', f_a, N(v(1800,fv=1.,d=.5,nu=20.,tabv=D_i,tab=D_i *-1)))
w.write('porta_fecha8.wav', f_a, N(v(1800,fv=1.,d=.5,nu=20.,tabv=D_i,tab=Q_i *-1)))


w.write('clique.wav', f_a, N(n.array([0]*100+[1]+[0]*10000)))
w.write('clique2.wav', f_a, N(adsr(v(fv=20,d=.2),S=-3.)))
w.write('clique3.wav', f_a, N(adsr(v(fv=20,d=.2,tab=Tr_i),S=-3.)))
w.write('clique4.wav', f_a, N(adsr(v(f=1000.,fv=20,d=.2,tab=Tr_i),S=-3.)))
w.write('clique5.wav', f_a, N(adsr(v(f=660.,fv=20,d=.2,tab=Tr_i),S=-3.)))


w.write('seleciona.wav', f_a, N(adsr(v(f=460.,fv=1.,d=.1,tab=Tr_i),S=-3.,R=10.)))
w.write('seleciona2.wav', f_a, N(adsr(v(f=460.,fv=10.,d=.1,tab=Tr_i),S=-3.,R=10.)))
w.write('cancela.wav', f_a, N(adsr(v(f=460.,fv=100.,d=.1,tab=Tr_i),S=-3.,R=10.)))
w.write('cancela2.wav', f_a, N(adsr(v(f=40.,fv=100.,d=.1,tab=Tr_i),S=-3.,R=10.)))


w.write('msgPos.wav', f_a, N(H((
           adsr(v(f=440.*(3/4.),fv=0.,nu=0.,d=.1,tab=Tr_i),S=-3.,R=10.),
           adsr(v(f=440.,fv=0.,nu=0.,d=.1,tab=Tr_i),S=-3.,R=10.),

))))

w.write('msgNeg.wav', f_a, N(H((
           adsr(v(f=440.*(3/4.),fv=0.,nu=0.,d=.1,tab=Tr_i),S=-3.,R=10.),
           adsr(v(f=440.,fv=0.,nu=0.,d=.1,tab=Tr_i),S=-3.,R=10.),
           adsr(v(f=440.*(3/4.),fv=0.,nu=0.,d=.1,tab=Tr_i),S=-3.,R=10.),
))))


w.write('msgPos2.wav', f_a, N(H((
           adsr(v(f=840.*(3/4.),fv=0.,nu=0.,d=.1,tab=Tr_i),S=-3.,R=10.),
           adsr(v(f=840.,fv=0.,nu=0.,d=.1,tab=Tr_i),S=-3.,R=10.),

))))


w.write('msgNeg2.wav', f_a, N(H((
           adsr(v(f=840.*(3/4.),fv=0.,nu=0.,d=.1,tab=Tr_i),S=-3.,R=10.),
           adsr(v(f=840.,fv=0.,nu=0.,d=.1,tab=Tr_i),S=-3.,R=10.),
           adsr(v(f=840.*(3/4.),fv=0.,nu=0.,d=.1,tab=Tr_i),S=-3.,R=10.),
))))

w.write('msgNeg3.wav', f_a, N(H((
           adsr(v(f=840.,fv=0.,nu=0.,d=.1,tab=Tr_i),S=-3.,R=10.),
           adsr(v(f=840.*(3/4.),fv=0.,nu=0.,d=.1,tab=Tr_i),S=-3.,R=10.),
))))


w.write('msgPos3.wav', f_a, N(H((
           adsr(v(f=840.*(3/4.),fv=0.,nu=0.,d=.1,tab=Tr_i),S=-3.,R=10.),
           adsr(v(f=840.,fv=0.,nu=0.,d=.1,tab=Tr_i),S=-3.,R=10.),
           adsr(v(f=840.*(2**(4./12)),fv=0.,nu=0.,d=.1,tab=Tr_i),S=-3.,R=10.),
           adsr(v(f=840.,fv=0.,nu=0.,d=.1,tab=Tr_i),S=-3.,R=10.),

))))

w.write('msgPos4.wav', f_a, N(H((
           adsr(v(f=840.*(3/4.),fv=0.,nu=0.,d=.1,tab=Tr_i),S=-3.,R=10.),
           adsr(v(f=840.,fv=0.,nu=0.,d=.1,tab=Tr_i),S=-3.,R=10.),
           adsr(v(f=840.*(2**(4./12)),fv=0.,nu=0.,d=.1,tab=Tr_i),S=-3.,R=10.),
           adsr(v(f=840.*(2**(7./12)),fv=0.,nu=0.,d=.1,tab=Tr_i),S=-3.,R=10.),
           adsr(v(f=840.,fv=0.,nu=0.,d=.1,tab=Tr_i),S=-3.,R=10.),

))))

w.write('msgNeg4.wav', f_a, N(H((
           adsr(v(f=840.,fv=0.,nu=0.,d=.1,tab=Tr_i),S=-3.,R=10.),
           adsr(v(f=840.*(2**(-6./12)),fv=0.,nu=0.,d=.1,tab=Tr_i),S=-3.,R=10.),
))))

w.write('perda.wav', f_a, N(H((
           adsr(v(f=840.,fv=0.,nu=0.,d=.1,              tab=D_i),S=-3.,R=10.),
           adsr(v(f=840.*(2**(-6./12)),fv=0.,nu=0.,d=.1,tab=D_i),S=-3.,R=10.),
))))

w.write('ganho.wav', f_a, N(H((
           adsr(v(f=840.*(2**(-7./12)),fv=0.,nu=0.,d=.1,tab=D_i),S=-3.,R=10.),
           adsr(v(f=840.,fv=0.,nu=0.,d=.1,              tab=D_i),S=-3.,R=10.),
))))

w.write('ganho2.wav', f_a, N(H((
           adsr(v(f=840.,fv=0.,nu=0.,d=.075,              tab=D_i),S=-3.,R=10.,A=5.,D=10.),
           adsr(v(f=840.,fv=0.,nu=0.,d=.025,              tab=D_i),S=-3.,R=10.,A=5.,D=10.),
           adsr(v(f=840.,fv=0.,nu=0.,d=.05,              tab=D_i),S=-3.,R=10.,A=5.,D=10.),
           adsr(v(f=840.,fv=0.,nu=0.,d=.05,              tab=D_i),S=-3.,R=5.,A=5.,D=10.),
))))

w.write('ganho3.wav', f_a, N(H((
           adsr(v(f=240.,fv=0.,nu=0.,d=.75,              tab=D_i),S=-9.,R=10.,A=5.,D=610.),
           adsr(v(f=240.*(2.**(-7/12.)),fv=0.,nu=0.,d=.25,              tab=D_i),S=-9.,R=10.,A=5.,D=210.),
           adsr(v(f=240.*(2.**(4/12.)),fv=0.,nu=0.,d=.5,              tab=D_i), S=-9.,R=10.,A=5., D=410.),
           adsr(v(f=240.*(2.**(0/12.)),fv=0.,nu=0.,d=.5,              tab=D_i), S=-9.,R=5.,A=5.,  D=410.),
))))

w.write('ganho4.wav', f_a, N(H((
           adsr(v(f=240.,fv=0.,nu=0.,d=.175,              tab=D_i),S=-9.,R=10.,A=5.,D=60.),
           adsr(v(f=240.*(2.**(-7/12.)),fv=0.,nu=0.,d=.25,              tab=D_i),S=-9.,R=10.,A=5.,D=210.),
           adsr(v(f=240.*(2.**(4/12.)),fv=0.,nu=0.,d=.105,              tab=D_i), S=-9.,R=10.,A=5., D=40.),
           adsr(v(f=240.*(2.**(0/12.)),fv=0.,nu=0.,d=.105,              tab=D_i), S=-9.,R=5.,A=5.,  D=40.),
))))

w.write('ganho5.wav', f_a, N(H((
           adsr(v(f=240.,fv=0.,nu=0.,d=.175,              tab=D_i),S=-9.,R=10.,A=5.,D=60.),
           adsr(v(f=240.*(2.**(7/12.)),fv=0.,nu=0.,d=.25,              tab=D_i),S=-9.,R=10.,A=5.,D=210.),
           adsr(v(f=240.*(2.**(4/12.)),fv=0.,nu=0.,d=.105,              tab=D_i), S=-9.,R=10.,A=5., D=40.),
           adsr(v(f=240.*(2.**(0/12.)),fv=0.,nu=0.,d=.105,              tab=D_i), S=-9.,R=5.,A=5.,  D=40.),
))))

w.write('ganho6.wav', f_a, N(H((
           adsr(v(f=240.,fv=0.,nu=0.,d=.175,              tab=D_i),S=-9.,R=10.,A=5.,D=60.),
           adsr(v(f=240.*(2.**(7/12.)),fv=0.,nu=0.,d=.25,              tab=D_i),S=-9.,R=10.,A=5.,D=210.),
           adsr(v(f=240.*(2.**(4/12.)),fv=0.,nu=0.,d=.105,              tab=D_i), S=-9.,R=10.,A=5., D=40.),
           adsr(v(f=240.*(2.**(12./12.)),fv=0.,nu=0.,d=.105,              tab=D_i), S=-9.,R=5.,A=5.,  D=40.),
))))


w.write('perda2.wav', f_a, N(H((
           adsr(v(f=240.,fv=0.,nu=0.,d=.25,              tab=D_i),S=-9.,R=10.,A=5.,D=60.)+
           adsr(v(f=240.*(2.**(6/12.)),fv=0.,nu=0.,d=.25,              tab=D_i),S=-9.,R=10.,A=5.,D=210.),
           adsr(v(f=240.*(2.**(4/12.)),fv=0.,nu=0.,d=.105,              tab=D_i), S=-9.,R=10.,A=5., D=40.)+
           adsr(v(f=240.*(2.**(3./12.)),fv=0.,nu=0.,d=.105,              tab=D_i), S=-9.,R=5.,A=5.,  D=40.),
))))


w.write('perda3.wav', f_a, N(H((
           adsr(v(f=240.,fv=0.,nu=0.,d=.25,              tab=D_i),S=-9.,R=10.,A=5.,D=60.)+
           adsr(v(f=240.*(2.**(6/12.)),fv=0.,nu=0.,d=.25,              tab=D_i),S=-9.,R=10.,A=5.,D=210.),
))))

w.write('perda4.wav', f_a, N(H((
           adsr(v(f=240.*(2.**(4/12.)),fv=0.,nu=0.,d=.105,              tab=D_i), S=-9.,R=10.,A=5., D=40.)+
           adsr(v(f=240.*(2.**(3./12.)),fv=0.,nu=0.,d=.105,              tab=D_i), S=-9.,R=5.,A=5.,  D=40.),

))))

w.write('perda5.wav', f_a, N(H((
           adsr(v(f=240.*(2.**(0/12.)),fv=0.,nu=0.,d=.105,              tab=D_i), S=-9.,R=10.,A=5., D=40.)+
           adsr(v(f=240.*(2.**(11./12.)),fv=0.,nu=0.,d=.105,              tab=D_i), S=-9.,R=5.,A=5.,  D=40.),

))))

w.write('ganhoX.wav', f_a, N(H((
           adsr(v(f=240.*(2.**(0/12.)),fv=0.,nu=0.,  d=.065,              tab=Q_i), S=-9.,R=10.,A=5., D=40.),
           adsr(v(f=240.*(2.**(11./12.)),fv=0.,nu=0.,d=.065,              tab=Q_i), S=-9.,R=5.,A=5.,  D=40.),
           adsr(v(f=240.*(2.**(0/12.)),fv=0.,nu=0.,  d=.065,              tab=Q_i), S=-9.,R=10.,A=5., D=40.),
           adsr(v(f=240.*(2.**(0/12.)),fv=0.,nu=0.,  d=.065,              tab=Q_i), S=-9.,R=10.,A=5., D=40.),
           adsr(v(f=240.*(2.**(0/12.)),fv=0.,nu=0.,  d=.065,              tab=Q_i), S=-9.,R=10.,A=5., D=40.),
           adsr(v(f=240.*(2.**(11./12.)),fv=0.,nu=0.,d=.065,              tab=Q_i), S=-9.,R=5.,A=5.,  D=40.),
           adsr(v(f=240.*(2.**(11./12.)),fv=0.,nu=0.,d=.065,              tab=D_i), S=-9.,R=5.,A=5.,  D=40.),
           adsr(v(f=240.*(2.**(11./12.)),fv=0.,nu=0.,d=.065,              tab=Q_i), S=-9.,R=5.,A=5.,  D=40.),

))))

w.write('ganhoX2.wav', f_a, N(H((
           adsr(v(f=240.*(2.**(0/12.)),fv=0.,nu=0.,  d=.065,              tab=Q_i), S=-9.,R=10.,A=5., D=40.),
           adsr(v(f=240.*(2.**(11./12.)),fv=0.,nu=0.,d=.065,              tab=Q_i), S=-9.,R=5.,A=5.,  D=40.),
           adsr(v(f=240.*(2.**(6/12.)),fv=0.,nu=0.,  d=.065,              tab=Q_i), S=-9.,R=10.,A=5., D=40.),
           adsr(v(f=240.*(2.**(0/12.)),fv=0.,nu=0.,  d=.065,              tab=Q_i), S=-9.,R=10.,A=5., D=40.),
           adsr(v(f=240.*(2.**(-1/12.)),fv=0.,nu=0.,  d=.065,              tab=Q_i), S=-9.,R=10.,A=5., D=40.),
           adsr(v(f=240.*(2.**(11./12.)),fv=0.,nu=0.,d=.065,              tab=Q_i), S=-9.,R=5.,A=5.,  D=40.),
           adsr(v(f=240.*(2.**(6./12.)),fv=0.,nu=0.,d=.065,              tab=D_i), S=-9.,R=5.,A=5.,  D=40.),
           adsr(v(f=240.*(2.**(0./12.)),fv=0.,nu=0.,d=.065,              tab=Q_i), S=-9.,R=5.,A=5.,  D=40.),

))))

w.write('ganhoX3.wav', f_a, N(H((
           adsr(v(f=240.*(2.**(0/12.)),fv=0.  ,nu=10.,  d=.065,              tab=Q_i), S=-9.,R=10.,A=5., D=40.),
           adsr(v(f=240.*(2.**(11./12.)),fv=0.,nu=10.,d=.065,              tab=Q_i), S=-9.,R=5.,A=5.,  D=40.),
           adsr(v(f=240.*(2.**(6/12.)),fv=0.  ,nu=10.,  d=.065,              tab=Q_i), S=-9.,R=10.,A=5., D=40.),
           adsr(v(f=240.*(2.**(0/12.)),fv=0.  ,nu=10.,  d=.065,              tab=Q_i), S=-9.,R=10.,A=5., D=40.),
           adsr(v(f=240.*(2.**(-1/12.)),fv=0. ,nu=10.,  d=.065,              tab=Q_i), S=-9.,R=10.,A=5., D=40.),
           adsr(v(f=240.*(2.**(11./12.)),fv=0.,nu=10.,d=.065,              tab=Q_i), S=-9.,R=5.,A=5.,  D=40.),
           adsr(v(f=240.*(2.**(6./12.)),fv=0. ,nu=10.,d=.065,              tab=D_i), S=-9.,R=5.,A=5.,  D=40.),
           adsr(v(f=240.*(2.**(0./12.)),fv=0. ,nu=10.,d=.065,              tab=Q_i), S=-9.,R=5.,A=5.,  D=40.),

))))

w.write('perdaX4.wav', f_a, N(H((
           adsr(v(f=240.*(2.**(0/12.)) , fv=100.  ,nu=10.,  d=.065,              tab=Q_i), S=-9.,R=10.,A=5., D=40.),
           adsr(v(f=240.*(2.**(11./12.)),fv=100.,nu=10.,d=.065,              tab=Q_i), S=-9.,R=5.,A=5.,  D=40.),
           adsr(v(f=240.*(2.**(6/12.)),  fv=100.  ,nu=10.,  d=.065,              tab=Q_i), S=-9.,R=10.,A=5., D=40.),
           adsr(v(f=240.*(2.**(0/12.)) , fv=100.  ,nu=10.,  d=.065,              tab=Q_i), S=-9.,R=10.,A=5., D=40.),
           adsr(v(f=240.*(2.**(-1/12.)), fv=100. ,nu=10.,  d=.065,              tab=Q_i), S=-9.,R=10.,A=5., D=40.),
           adsr(v(f=240.*(2.**(11./12.)),fv=100.,nu=10.,d=.065,              tab=Q_i), S=-9.,R=5.,A=5.,  D=40.),
           adsr(v(f=240.*(2.**(6./12.)), fv=100. ,nu=10.,d=.065,              tab=D_i), S=-9.,R=5.,A=5.,  D=40.),
           adsr(v(f=240.*(2.**(0./12.)), fv=100. ,nu=10.,d=.065,              tab=Q_i), S=-9.,R=5.,A=5.,  D=40.),

))))

w.write('perdaX5.wav', f_a, N(H((
           adsr(v(f=240.*(2.**(0/12.)) , fv=200.  ,nu=10.,  d=.065,              tab=Q_i), S=-9.,R=10.,A=5., D=40.),
           adsr(v(f=240.*(2.**(11./12.)),fv=200.,nu=10.,d=.065,              tab=Q_i), S=-9.,R=5.,A=5.,  D=40.),
           adsr(v(f=240.*(2.**(6/12.)),  fv=200.  ,nu=10.,  d=.065,              tab=Q_i), S=-9.,R=10.,A=5., D=40.),
           adsr(v(f=240.*(2.**(0/12.)) , fv=200.  ,nu=10.,  d=.065,              tab=Q_i), S=-9.,R=10.,A=5., D=40.),
           adsr(v(f=240.*(2.**(-1/12.)), fv=200. ,nu=10.,  d=.065,              tab=Q_i), S=-9.,R=10.,A=5., D=40.),
           adsr(v(f=240.*(2.**(11./12.)),fv=200.,nu=10.,d=.065,              tab=Q_i), S=-9.,R=5.,A=5.,  D=40.),
           adsr(v(f=240.*(2.**(6./12.)), fv=200. ,nu=10.,d=.065,              tab=D_i), S=-9.,R=5.,A=5.,  D=40.),
           adsr(v(f=240.*(2.**(0./12.)), fv=200. ,nu=10.,d=.065,              tab=Q_i), S=-9.,R=5.,A=5.,  D=40.),

))))

w.write('videogame.wav', f_a, N(H((
           adsr(v(f=240.*(2.**(0/12.)),fv=0.,nu=0.,  d=.065,              tab=Tr_i), S=-9.,R=10.,A=5., D=40.),
           adsr(v(f=240.*(2.**(4./12.)),fv=0.,nu=0.,d=.065,              tab=Tr_i), S=-9.,R=5.,A=5.,  D=40.),
           adsr(v(f=240.*(2.**(7/12.)),fv=0.,nu=0.,  d=.065,              tab=Tr_i), S=-9.,R=10.,A=5., D=40.),
           adsr(v(f=240.*(2.**(7/12.)),fv=0.,nu=0.,  d=.065,              tab=Tr_i), S=-9.,R=10.,A=5., D=40.),
           adsr(v(f=240.*(2.**(4/12.)),fv=0.,nu=0.,  d=.065,             tab=Tr_i), S=-9.,R=10.,A=5., D=40.),
           adsr(v(f=240.*(2.**(0./12.)),fv=0.,nu=0.,d=.065,              tab=Tr_i), S=-9.,R=5.,A=5.,  D=40.),
           adsr(v(f=240.*(2.**(7./12.)),fv=0.,nu=0.,d=.065,               tab=Tr_i), S=-9.,R=5.,A=5.,  D=40.),
           adsr(v(f=240.*(2.**(0./12.)),fv=0.,nu=0.,d=.065,               tab=Tr_i), S=-9.,R=5.,A=5.,  D=40.),

           adsr(v(f=240.*(2.**(0/12.)),fv=0.,nu=0.,  d=.65,             tab=Tr_i), S=-9.,R=10.,A=5., D=40.),
           adsr(v(f=240.*(2.**(11./12.)),fv=0.,nu=0.,d=.065,            tab=Tr_i), S=-9.,R=5.,A=5.,  D=40.),
           adsr(v(f=240.*(2.**(12/12.)),fv=0.,nu=0.,  d=.065,           tab=Tr_i), S=-9.,R=10.,A=5., D=40.),
           adsr(v(f=240.*(2.**(11/12.)),fv=0.,nu=0.,  d=.65,            tab=Tr_i), S=-9.,R=10.,A=5., D=40.),
           adsr(v(f=240.*(2.**(12/12.)),fv=0.,nu=0.,  d=.065,           tab=Tr_i), S=-9.,R=10.,A=5., D=40.),
           adsr(v(f=240.*(2.**(0./12.)),fv=0.,nu=0.,d=.65,              tab=Tr_i), S=-9.,R=5.,A=5.,  D=40.),
           adsr(v(f=240.*(2.**(6./12.)),fv=0.,nu=0.,d=.065,             tab=Tr_i), S=-9.,R=5.,A=5.,  D=40.),
           adsr(v(f=240.*(2.**(-1./12.)),fv=0.,nu=0.,d=.65,             tab=Tr_i), S=-9.,R=5.,A=5.,  D=40.),

           adsr(v(f=240.*(2.**(0/12.)),fv=0.,nu=0.,  d=.065,              tab=Tr_i), S=-9.,R=10.,A=5., D=40.),
           adsr(v(f=240.*(2.**(11./12.)),fv=0.,nu=0.,d=.065,              tab=Q_i), S=-9.,R=5.,A=5.,  D=40.),
           adsr(v(f=240.*(2.**(12/12.)),fv=0.,nu=0.,  d=.5,              tab=Tr_i), S=-9.,R=10.,A=5., D=40.),
           adsr(v(f=240.*(2.**(14/12.)),fv=0.,nu=0.,  d=.065,              tab=Q_i), S=-9.,R=10.,A=5., D=40.),
           adsr(v(f=240.*(2.**(16/12.)),fv=0.,nu=0.,  d=.5,              tab=Tr_i), S=-9.,R=10.,A=5., D=40.),
           adsr(v(f=240.*(2.**(17./12.)),fv=0.,nu=0.,d=.065,              tab=Q_i), S=-9.,R=5.,A=5.,  D=40.),
           adsr(v(f=240.*(2.**(19./12.)),fv=0.,nu=0.,d=.65,              tab=Tr_i), S=-9.,R=5.,A=5.,  D=40.),
           adsr(v(f=240.*(2.**(12./12.)),fv=0.,nu=0.,d=.065,              tab=Tr_i), S=-9.,R=5.,A=5.,  D=40.),
))))


w.write('videogame2.wav', f_a, N(H((
           adsr(v(f=240.*(2.**(0/12.)),fv=0.,nu=0.,  d=.065*4,              tab=Tr_i), S=-9.,R=10.,A=5., D=40.),
           adsr(v(f=240.*(2.**(4./12.)),fv=0.,nu=0.,d=.0652*2,              tab=Tr_i), S=-9.,R=5.,A=5.,  D=40.),
           adsr(v(f=240.*(2.**(7/12.)),fv=0.,nu=0.,  d=.065,              tab=Tr_i), S=-9.,R=10.,A=5., D=40.),
           adsr(v(f=240.*(2.**(7/12.)),fv=0.,nu=0.,  d=.065,              tab=Tr_i), S=-9.,R=10.,A=5., D=40.),
           adsr(v(f=240.*(2.**(4/12.)),fv=0.,nu=0.,  d=.065,             tab=Tr_i), S=-9.,R=10.,A=5., D=40.),
           adsr(v(f=240.*(2.**(0./12.)),fv=0.,nu=0.,d=.065,              tab=Tr_i), S=-9.,R=5.,A=5.,  D=40.),
           adsr(v(f=240.*(2.**(7./12.)),fv=0.,nu=0.,d=.065,               tab=Tr_i), S=-9.,R=5.,A=5.,  D=40.),
           adsr(v(f=240.*(2.**(0./12.)),fv=0.,nu=0.,d=.065,               tab=Tr_i), S=-9.,R=5.,A=5.,  D=40.),

           adsr(v(f=240.*(2.**(0/12.)),fv=0.,nu=0.,  d=.065,             tab=Tr_i), S=-9.,R=10.,A=5., D=40.),
           adsr(v(f=240.*(2.**(11./12.)),fv=0.,nu=0.,d=.065,            tab=Tr_i), S=-9.,R=5.,A=5.,  D=40.),
           adsr(v(f=240.*(2.**(12/12.)),fv=0.,nu=0.,  d=.065,           tab=Tr_i), S=-9.,R=10.,A=5., D=40.),
           adsr(v(f=240.*(2.**(11/12.)),fv=0.,nu=0.,  d=.065*4,            tab=Tr_i), S=-9.,R=10.,A=5., D=40.),
           adsr(v(f=240.*(2.**(12/12.)),fv=0.,nu=0.,  d=.065,           tab=Tr_i), S=-9.,R=10.,A=5., D=40.),
           adsr(v(f=240.*(2.**(0./12.)),fv=0.,nu=0.,d=.065,              tab=Tr_i), S=-9.,R=5.,A=5.,  D=40.),
           adsr(v(f=240.*(2.**(6./12.)),fv=0.,nu=0.,d=.065,             tab=Tr_i), S=-9.,R=5.,A=5.,  D=40.),
           adsr(v(f=240.*(2.**(-1./12.)),fv=0.,nu=0.,d=.065*2,             tab=Tr_i), S=-9.,R=5.,A=5.,  D=40.),

           adsr(v(f=240.*(2.**(0/12.)),fv=0.,nu=0.,  d=.065,              tab=Tr_i), S=-9.,R=10.,A=5., D=40.),
           adsr(v(f=240.*(2.**(11./12.)),fv=0.,nu=0.,d=.065,              tab=Tr_i), S=-9.,R=5.,A=5.,  D=40.),
           adsr(v(f=240.*(2.**(12/12.)),fv=0.,nu=0.,  d=.065,              tab=Tr_i), S=-9.,R=10.,A=5., D=40.),
           adsr(v(f=240.*(2.**(14/12.)),fv=0.,nu=0.,  d=.065*4,              tab=Tr_i), S=-9.,R=10.,A=5., D=40.),
           adsr(v(f=240.*(2.**(16/12.)),fv=0.,nu=0.,  d=.065,              tab=Tr_i), S=-9.,R=10.,A=5., D=40.),
           adsr(v(f=240.*(2.**(17./12.)),fv=0.,nu=0.,d=.065,              tab=Tr_i), S=-9.,R=5.,A=5.,  D=40.),
           adsr(v(f=240.*(2.**(19./12.)),fv=0.,nu=0.,d=.065*2,              tab=Tr_i), S=-9.,R=5.,A=5.,  D=40.),
           adsr(v(f=240.*(2.**(12./12.)),fv=0.,nu=0.,d=.065,              tab=Tr_i), S=-9.,R=5.,A=5.,  D=40.),
))))


w.write('videogame3.wav', f_a, N(H((
           adsr(v(f=240.*(2.**(0/12.)),fv=0.,nu=0.,  d=.065*4,              tab=Tr_i), S=-9.,R=10.,A=5., D=40.),
           adsr(v(f=240.*(2.**(4./12.)),fv=0.,nu=0.,d=.0652*2,              tab=Tr_i), S=-9.,R=5.,A=5.,  D=40.),
           adsr(v(f=240.*(2.**(7/12.)),fv=0.,nu=0.,  d=.065,              tab=Tr_i), S=-9.,R=10.,A=5., D=40.),
           adsr(v(f=240.*(2.**(7/12.)),fv=0.,nu=0.,  d=.065,              tab=Tr_i), S=-9.,R=10.,A=5., D=40.),
           adsr(v(f=240.*(2.**(4/12.)),fv=0.,nu=0.,  d=.065,             tab=Tr_i), S=-9.,R=10.,A=5., D=40.),
           adsr(v(f=240.*(2.**(0./12.)),fv=0.,nu=0.,d=.065*3,              tab=Tr_i), S=-9.,R=5.,A=5.,  D=40.),
           adsr(v(f=240.*(2.**(7./12.)),fv=0.,nu=0.,d=.065,               tab=Tr_i), S=-9.,R=5.,A=5.,  D=40.),
           adsr(v(f=240.*(2.**(0./12.)),fv=0.,nu=0.,d=.065,               tab=Tr_i), S=-9.,R=5.,A=5.,  D=40.),

           adsr(v(f=240.*(2.**(0/12.)),fv=0.,nu=0.,  d=.065,             tab=Tr_i), S=-9.,R=10.,A=5., D=40.),
           adsr(v(f=240.*(2.**(11./12.)),fv=0.,nu=0.,d=.065*3,            tab=Tr_i), S=-9.,R=5.,A=5.,  D=40.),
           adsr(v(f=240.*(2.**(12/12.)),fv=0.,nu=0.,  d=.065,           tab=Tr_i), S=-9.,R=10.,A=5., D=40.),
           adsr(v(f=240.*(2.**(11/12.)),fv=0.,nu=0.,  d=.065*4,            tab=Tr_i), S=-9.,R=10.,A=5., D=40.),
           adsr(v(f=240.*(2.**(12/12.)),fv=0.,nu=0.,  d=.065,           tab=Tr_i), S=-9.,R=10.,A=5., D=40.),
           adsr(v(f=240.*(2.**(0./12.)),fv=0.,nu=0.,d=.065,              tab=Tr_i), S=-9.,R=5.,A=5.,  D=40.),
           adsr(v(f=240.*(2.**(6./12.)),fv=0.,nu=0.,d=.065,             tab=Tr_i), S=-9.,R=5.,A=5.,  D=40.),
           adsr(v(f=240.*(2.**(-1./12.)),fv=0.,nu=0.,d=.065*2,             tab=Tr_i), S=-9.,R=5.,A=5.,  D=40.),

           adsr(v(f=240.*(2.**(0/12.)),fv=0.,nu=0.,  d=.065*3,              tab=Tr_i), S=-9.,R=10.,A=5., D=40.),
           adsr(v(f=240.*(2.**(11./12.)),fv=0.,nu=0.,d=.065,              tab=Tr_i), S=-9.,R=5.,A=5.,  D=40.),
           adsr(v(f=240.*(2.**(12/12.)),fv=0.,nu=0.,  d=.065,              tab=Tr_i), S=-9.,R=10.,A=5., D=40.),
           adsr(v(f=240.*(2.**(14/12.)),fv=0.,nu=0.,  d=.065*4,              tab=Tr_i), S=-9.,R=10.,A=5., D=40.),
           adsr(v(f=240.*(2.**(16/12.)),fv=0.,nu=0.,  d=.065,              tab=Tr_i), S=-9.,R=10.,A=5., D=40.),
           adsr(v(f=240.*(2.**(17./12.)),fv=0.,nu=0.,d=.065,              tab=Tr_i), S=-9.,R=5.,A=5.,  D=40.),
           adsr(v(f=240.*(2.**(19./12.)),fv=0.,nu=0.,d=.065*2,              tab=Tr_i), S=-9.,R=5.,A=5.,  D=40.),
           adsr(v(f=240.*(2.**(12./12.)),fv=0.,nu=0.,d=.065,              tab=Tr_i), S=-9.,R=5.,A=5.,  D=40.),
))))

w.write('videogame4.wav', f_a, N(H((
           adsr(v(f=240.*(2.**(0/12.)),fv=0.,nu=0.,  d=.065*4,              tab=Tr_i), S=-9.,R=10.,A=5., D=40.),
           adsr(v(f=240.*(2.**(4./12.)),fv=0.,nu=0.,d=.0652*2,              tab=Tr_i), S=-9.,R=5.,A=5.,  D=40.),
           adsr(v(f=240.*(2.**(7/12.)),fv=0.,nu=0.,  d=.065,              tab=Tr_i), S=-9.,R=10.,A=5., D=40.),
           adsr(v(f=240.*(2.**(7/12.)),fv=0.,nu=0.,  d=.065,              tab=Tr_i), S=-9.,R=10.,A=5., D=40.),
           adsr(v(f=240.*(2.**(4/12.)),fv=0.,nu=0.,  d=.065,             tab=Tr_i), S=-9.,R=10.,A=5., D=40.),
           adsr(v(f=240.*(2.**(0./12.)),fv=0.,nu=0.,d=.065*3,              tab=Tr_i), S=-9.,R=5.,A=5.,  D=40.),
           adsr(v(f=240.*(2.**(7./12.)),fv=0.,nu=0.,d=.065,               tab=Tr_i), S=-9.,R=5.,A=5.,  D=40.),
           adsr(v(f=240.*(2.**(0./12.)),fv=0.,nu=0.,d=.065,               tab=Tr_i), S=-9.,R=5.,A=5.,  D=40.),

           adsr(v(f=240.*(2.**(0/12.)),fv=0.,nu=0.,  d=.065*5,             tab=Tr_i), S=-9.,R=10.,A=5., D=40.),
           adsr(v(f=240.*(2.**(11./12.)),fv=0.,nu=0.,d=.065*3,            tab=Tr_i), S=-9.,R=5.,A=5.,  D=40.),
           adsr(v(f=240.*(2.**(12/12.)),fv=0.,nu=0.,  d=.065,           tab=Tr_i), S=-9.,R=10.,A=5., D=40.),
           adsr(v(f=240.*(2.**(11/12.)),fv=0.,nu=0.,  d=.065*4,            tab=Tr_i), S=-9.,R=10.,A=5., D=40.),
           adsr(v(f=240.*(2.**(12/12.)),fv=0.,nu=0.,  d=.065,           tab=Tr_i), S=-9.,R=10.,A=5., D=40.),
           adsr(v(f=240.*(2.**(0./12.)),fv=0.,nu=0.,d=.065,              tab=Tr_i), S=-9.,R=5.,A=5.,  D=40.),
           adsr(v(f=240.*(2.**(2./12.)),fv=0.,nu=0.,d=.065,             tab=Tr_i), S=-9.,R=5.,A=5.,  D=40.),
           adsr(v(f=240.*(2.**(-1./12.)),fv=0.,nu=0.,d=.065*2,             tab=Tr_i), S=-9.,R=5.,A=5.,  D=40.),

           adsr(v(f=240.*(2.**(0/12.)),fv=0.,nu=0.,  d=.065*3,              tab=Tr_i), S=-9.,R=10.,A=5., D=40.),
           adsr(v(f=240.*(2.**(11./12.)),fv=0.,nu=0.,d=.065,              tab=Tr_i), S=-9.,R=5.,A=5.,  D=40.),
           adsr(v(f=240.*(2.**(12/12.)),fv=0.,nu=0.,  d=.065*4,              tab=Tr_i), S=-9.,R=10.,A=5., D=40.),
           adsr(v(f=240.*(2.**(14/12.)),fv=0.,nu=0.,  d=.065*4,              tab=Tr_i), S=-9.,R=10.,A=5., D=40.),
           adsr(v(f=240.*(2.**(16/12.)),fv=0.,nu=0.,  d=.065,              tab=Tr_i), S=-9.,R=10.,A=5., D=40.),
           adsr(v(f=240.*(2.**(17./12.)),fv=0.,nu=0.,d=.065,              tab=Tr_i), S=-9.,R=5.,A=5.,  D=40.),
           adsr(v(f=240.*(2.**(19./12.)),fv=0.,nu=0.,d=.065*2,              tab=Tr_i), S=-9.,R=5.,A=5.,  D=40.),
           adsr(v(f=240.*(2.**(12./12.)),fv=0.,nu=0.,d=.065,              tab=Tr_i), S=-9.,R=5.,A=5.,  D=40.),
))))



# abre todos os gritoFala*
# passa por um passa bandas que soh passa uns medios
# salva como tv_gritoFala*

#
#c = n.zeros(len(coefs))
#c[1000:10000] = n.exp(1j*n.random.uniform(0, 2*n.pi, 9000))
#
## real par, imaginaria impar
#c[Lambda/2+1:] = n.real(c[1:Lambda/2])[::-1] - 1j * \
#    n.imag(c[1:Lambda/2])[::-1]
#
#resp_imp= n.fft.ifft(c)
#resp_imp_= n.real(resp_imp)
#import os
#
#ll=os.listdir(".")
#ll=[lll for lll in ll if "gritoFala" in lll]
#for i in ll:
#    print i
#    foo=n.convolve(w.read("%s"%(i,))[1],resp_imp)
#    w.write('tv_%s'%(i,),  f_a, N(foo))
#    print i
#









