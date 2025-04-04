import numpy as np
import pandas as pd
import sympy as sp

class Funciones:

  def __init__(self) -> None:
     self.x = sp.Symbol('x')

  #Biseccion
  def Biseccion(self, f, a, b, tol):
    if f(a)*f(b) > 0:
      print("la funcion no cumple el teorema en el intervalo dado")
      return
    else:
      print("Buscando.....")
      i=0
      while (abs(b-a)>tol):
        p = (a+b)/2
        i+=1
        if(f(a)*f(p)>0):
          a=p
        else:
          b=p
      return p, i, abs(b-a), "Biseccion"


  #Falsa posicion
  def Falsa_posicion(self, f, a, b, tol):
    if f(a)*f(b) > 0:
      print("la funcion no cumple el teorema en el intervalo dado")
      return
    else:
      print("Buscando.....")
      i=0
      error = abs(b - a)
      p = b-(f(b)*(a-b))/(f(a)-f(b))
      while error > tol:
        p = b-(f(b)*(a-b))/(f(a)-f(b))
        i+=1
        if(f(a)*f(p)>0):
          a=p
        else:
          b=p
        error = abs(b - a)
      return p, i, abs(b-a), "Falsa Posicion"


  #Newton
  def newton(self,f,x0, tol):
    f_expr = sp.sympify(f(self.x))
    f_prima = sp.diff(f_expr, self.x)  #Primera derivada
    print("Derivada:", f_prima)

    x1 = x0 - (f_expr.subs(self.x, x0) / f_prima.subs(self.x, x0))
    i=0
    error=abs(x1 - x0)
    while error > tol:
      x0 = x1
      x1 = x0 - sp.N((f_expr.subs(self.x, x0) / f_prima.subs(self.x, x0)))
      i+=1
      error = abs(x1 - x0)
    return round(x1, 6), i, round(error, 6), "Newton"


  #Secante
  def secante(self, f, x0, x1, tol):
    error = 1
    i=0
    while error > tol:
      x2= x1 - f(x1) * (x0 - x1) / (f(x0) - f(x1))
      i+=1
      error = abs(x2 - x1)
      x0 = x1
      x1 = x2
    return x2, i, error, "Secante"

  #Newton_detalle
  def newton_detalle(self,f,x0, tol):
    f_expr = sp.sympify(f(self.x))
    f_prima = sp.diff(f_expr, self.x)  #Primera derivada
    print("Derivada:", f_prima)

    x1 = x0 - (f_expr.subs(self.x, x0) / f_prima.subs(self.x, x0))
    i=0
    error=abs(x1 - x0)
    info = []
    while error > tol:
      x0 = x1
      x1 = x0 - sp.N((f_expr.subs(self.x, x0) / f_prima.subs(self.x, x0)))
      i+=1
      error = abs(x1 - x0)
      info.append([i, round(x0, 6), round(x1, 6), round(error, 6)])
    data = pd.DataFrame(info, columns=["Iteracion", "x0", "x1", "Error"])
    return round(x1, 6), data


