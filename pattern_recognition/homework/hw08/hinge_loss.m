function resid = hinge_loss(x, b, A, lambda)

    resid = sum(1 - b * A.T * x) + lambda * norm(x);

end
