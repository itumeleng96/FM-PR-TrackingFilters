from numpy import dot


def kf_predict(X,P,A,Q,U):

    """
    X: The mean state estimate of the previous step (1k−).
    P: The state covariance of previous step (1k−).
    A: The transition n    n×matrix.
    Q: The process noise covariance matrix.
    B: The input effect matrix.
    U: The control input
    """
    x=dot(A,X) +dot(B,U)
    P=dot(A,dot(P,A.T))+Q

    return(X,P)

def kf_update(X,P,Y,H,R):
    IM = dot(H,X)
    IS = R + dot(H,dot(P,H.T))
    K = dot(P,dot(H.T,inv(IS)))
    X = X + dot(K,(Y-IM))
    P = P - dot(K,dot(Y,IM,IS))
    LH = gauss_pdf(Y,IM,IS)

    return (X,P,K,IM,IS,LH)

def gauss_pdf(X,M,S):
    
