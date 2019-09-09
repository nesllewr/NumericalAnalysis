# <b> _Background_ </b>

- 머신 엡실론(machine epsilon) : 반올림 오차의 상한값으로 1과 1과 구별될 수 있는 1보다 큰 다음 수 사이의 간격  
-> 1과 1 + eps 사이에는 어떤 수도 나타낼 수 없다. 
- 한정된 자리수로 인해 컴퓨터 연산에서 오차가 발생한다.  



![na_hw1정의](https://user-images.githubusercontent.com/41321080/64501729-8ca9ae00-d2fd-11e9-93ab-ac9c929cfbbf.PNG)

# <b> _CODE_ </b>  
##            **main** 
```c
int main() {

	//1번째 방법
	printf("Method 1 (machar)\n");
	int ibeta, it, irnd, ngrd, machep, negep, iexp, minexp, maxexp;

	float eps, epsneg, xmin, xmax;
	machar(&ibeta, &it, &irnd, &ngrd, &machep, &negep, &iexp, &minexp, &maxexp, &eps, &epsneg, &xmin, &xmax);
	printf("float:  %g (%.23f)\n", eps, eps);

	double epsD, epsnegD, xminD, xmaxD;
	macharD(&ibeta, &it, &irnd, &ngrd, &machep, &negep, &iexp, &minexp, &maxexp, &epsD, &epsnegD, &xminD, &xmaxD);
	printf("double: %g (%.52f)\n", epsD, epsD);

	printf("\n");

	//2번째 방법
	printf("Method 2 (get_eps)\n");
	printf("double : %g (%.52f)\n", double_get_eps(), double_get_eps());
	printf("float : %g (%.23f)\n",float_get_eps(), float_get_eps());

	printf("\n");
	return 0;
}
```

- 결과값의 출력은 %g (%(최대자리수f))로 실수형 자동출력과 실수형 둘 다 보여준다.  
* method 1 : machar()함수를 사용해 epsilon을 구한다
* method 2  : float와 double 각각 1 + 2 ^(-n) = 1을 만독하는 최소 n를 구하는 get_eps()의 리턴값을 보여준다.  
 

---

## **function**

```c
void machar(int *ibeta, int *it, int *irnd, int *ngrd, int *machep, int *negep,
	int *iexp, int *minexp, int *maxexp, float *eps, float *epsneg,
	float *xmin, float *xmax);
void macharD(int *ibeta, int *it, int *irnd, int *ngrd, int *machep, int *negep,
	int *iexp, int *minexp, int *maxexp, double *eps, double *epsneg,
	double *xmin, double *xmax);

```
-macharD()는 machar()에서 float으로 선언된 부분을 double형으로 바꾼 함수


```c
double double_get_eps() {

	double eps = 1;
	while ((1 + eps / 2) > 1)
		eps /= 2;

	return eps;
}
 ```
```c
 float float_get_eps() {

	float eps = 1.0f;

	while (1.0f + eps / 2.0f > 1.0f)
		eps /= 2.0f;

	return eps;
}
 ```


-구하려는 자료형에 맞게 eps을 선언한다.    
-1 + 2^(-n) = 1 에서 2^(-n)부분을 eps/2로 본다. (0보다 큰 n을 제곱하는 것은 2로 n번 나누는 것과 동일)  
-while문을 돌면서 조건문 좌변의 값이 1보다 크지 않으면 더 이상 차이를 표현할 수 없어진다.  
-> 더 이상 차이를 표현할 수 없는 값이므로 각 자료형의 epsilon으로 값을 리턴한다.  


---

# <b>_Execution & Result_</b>
  
gcc (Ubuntu 5.4.0-6ubuntu1~16.04.10) 5.4.0 20160609  

![na_hw1결과2](https://user-images.githubusercontent.com/41321080/64501569-da71e680-d2fc-11e9-9c72-5b7bf8003310.PNG)
