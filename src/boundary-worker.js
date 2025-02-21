// boundary-worker.js
importScripts("https://cdnjs.cloudflare.com/ajax/libs/mathjs/10.6.4/math.js");

onmessage = function (event) {
    const { boundaryGaussians, intensityThreshold, removeCenter, removeRadius } = event.data;

    // Here, you'd put the logic to process each gaussian, e.g.:
    //  1) approximateGaussianOutsideCube(...)
    //  2) gather newly created gaussians
    //  3) return them

    const newEllipsoids = {
        positions: [],
        cov3Da: [],
        cov3Db: [],
        colors: [],
        opacities: [],
    };


    for (let g_indx = 0; g_indx < boundaryGaussians.length; g_indx++) {

        const subdivided = approximateGaussianOutsideCube(
            g_indx,
            boundaryGaussians,
            intensityThreshold,
            removeCenter,
            removeRadius
        );

        // Collect them
        newEllipsoids.positions.push(...subdivided.positions);
        newEllipsoids.cov3Da.push(...subdivided.cov3Da);
        newEllipsoids.cov3Db.push(...subdivided.cov3Db);
        newEllipsoids.colors.push(...subdivided.colors);
        newEllipsoids.opacities.push(subdivided.opacites);

    }

    // Done processing? Send results back:
    postMessage({ newEllipsoids });
};

// Option A: Define a global error handler:
self.onerror = function (message, source, lineno, colno, error) {
    // Forward the error details back to the main thread
    self.postMessage({
        type: "worker-error",
        message,
        source,
        lineno,
        colno,
        stack: error?.stack,
    });

    // Return true to prevent default logging (optional)
    return true;
};

function approximateGaussianOutsideCube(g_indx, boundaryGaussians, intensityThreshold, removeCenter, removeRadius) {
    // console.log("approximateGaussianOutsideCube");
    const resultEllipsoids = {
        positions: [],
        cov3Da: [],
        cov3Db: [],
        colors: [],
        opacities: [],
    };




    // Define the six planes of the cube
    const planes = [
        {
            normal: [0, 0, -1], point: [
                removeCenter[0],
                removeCenter[1],
                removeCenter[2]]
        }, // Bottom face
        {
            normal: [0, 0, 1], point: [
                removeCenter[0],
                removeCenter[1],
                removeCenter[2] + removeRadius * 2]
        }, // Top face
        {
            normal: [-1, 0, 0], point: [
                removeCenter[0] - removeRadius,
                removeCenter[1],
                removeCenter[2] + removeRadius]
        }, // Left face
        {
            normal: [1, 0, 0], point: [
                removeCenter[0] + removeRadius,
                removeCenter[1],
                removeCenter[2] + removeRadius]
        }, // Right face
        {
            normal: [0, -1, 0], point: [
                removeCenter[0],
                removeCenter[1] - removeRadius,
                removeCenter[2] + removeRadius]
        }, // Front face
        {
            normal: [0, 1, 0], point: [
                removeCenter[0],
                removeCenter[1] + removeRadius,
                removeCenter[2] + removeRadius]
        }  // Back face
    ];

    // Loop through each plane and call the helper function
    for (const { normal, point } of planes) {
        // approximate using 3 balls
        const smallerEllipsoids = approximateGaussianOutsideHalfSpace(g_indx, boundaryGaussians, intensityThreshold, point, normal, 3);


        smallerEllipsoids.forEach(ellipsoid => {
            resultEllipsoids.positions.push(...ellipsoid.position.toArray());
            resultEllipsoids.cov3Da.push(...ellipsoid.cov3Da);
            resultEllipsoids.cov3Db.push(...ellipsoid.cov3Db);
            resultEllipsoids.colors.push(...ellipsoid.color);
            resultEllipsoids.opacities.push(ellipsoid.opacity);
        });
    }

    return resultEllipsoids;
}


function approximateGaussianOutsideHalfSpace(g_indx, boundaryGaussians, intensityThreshold, planeCenter, planeNormal, numSmallerBalls) {
    // read data

    gaussian = boundaryGaussians[g_indx];
    const gPos = Array.from(gaussian.position);
    const [a, b, c] = gaussian.cov3Da;
    const [d, e, f] = gaussian.cov3Db;

    const gCov = math.matrix([
        [a, b, c],
        [b, d, e],
        [c, e, f],
    ]);

    const color = gaussian.color;
    const opacity = gaussian.opacity;

    // Eigen decomposition of Sigma to get U and Lambda
    const { values: eigenvalues, vectors: eigenvectors } = math.eigs(gCov);

    // Construct the affine transformation matrix A = U^T * Lambda^(-1/2) * U
    const LambdaInvSqrt = math.diag(eigenvalues.map(v => Math.sqrt(1 / v))); // Lambda^(-1/2)
    const U = eigenvectors;
    const A = math.multiply(math.multiply(math.transpose(U), LambdaInvSqrt), U); // A = U^T * Lambda^(-1/2) * U

    // Apply transformations
    const transformedPosition = math.multiply(A, gPos); // A * mu
    const transformedPlanePoint = math.multiply(A, planeCenter); // Transform removeCenter
    let transformedPlaneNormal = math.multiply(math.inv(A), planeNormal);
    transformedPlaneNormal = math.divide(transformedPlaneNormal, math.norm(transformedPlaneNormal)); // Normalize the normal

    // Calculate the determinant of the covariance matrix to determine scaling factor
    const scalingFactor = Math.sqrt((2 * Math.PI) ** 3 * Math.abs(math.det(gCov))) * intensityThreshold;
    const C = -2 * Math.log(scalingFactor);

    // Define the radius R of the transformed Gaussian (which is now a ball)
    const R = Math.sqrt(C); // Effective radius based on contour level

    // Use the utility function to get the smaller balls
    const ball = {
        center: transformedPosition,
        radius: R
    }
    let smallerBalls = []
    // console.log("Computing smaller balls...");
    if (numSmallerBalls === 2) {
        smallerBalls.push(...approximateCutBall2(ball, transformedPlanePoint, transformedPlaneNormal));
    }
    else if (numSmallerBalls === 3) {
        smallerBalls.push(...approximateCutBall3(ball, transformedPlanePoint, transformedPlaneNormal));
    }

    // Transform the smaller balls back to the original space
    const resultEllipsoids = smallerBalls.map(ball => {
        const newPos = math.multiply(math.inv(A), ball.center); // Transform back
        const newCov = math.multiply(gCov, ball.radius ** 2 / C)._data;
        const newCov3Da = [newCov[0][0], newCov[0][1], newCov[0][2]]; // Upper triangular part, row 0
        const newCov3Db = [newCov[1][1], newCov[1][2], newCov[2][2]]; // Upper triangular part, row 1
        return {
            position: newPos,
            cov3Da: newCov3Da,
            cov3Db: newCov3Db,
            color: color,
            opacity: opacity,
        };
    });

    return resultEllipsoids;
}

function approximateCutBall2(originalBall, planePoint, planeNormal) {
    // Calculate the distance h from the ball center to the plane
    const h = math.dot(math.subtract(originalBall.center, planePoint), planeNormal);

    // Check if the ball intersects the plane
    if (h >= originalBall.radius) {
        return [originalBall]; // Ball is entirely above the plane
    }

    // To find an offset vector parallel to the plane, we can pick one vector in the plane
    let planeTangent = math.cross(planeNormal, [1, 0, 0]); // First tangent vector

    // Ensure we have a valid tangent, recalculate if the first cross product resulted in a zero vector
    if (math.norm(planeTangent) === 0) {
        planeTangent = math.cross(planeNormal, [0, 1, 0]);
    }

    let r, centerPoint;

    if (h >= originalBall.radius / 2) {
        r = originalBall.radius / 2;
        centerPoint = originalBall.center;
    }
    else {
        // Solve for r using the quadratic equation: r^2 + (2R - 2h)r + (h^2 - R^2) = 0
        const a = 1; // Coefficient of r^2
        const b = 2 * originalBall.radius - 2 * h; // Coefficient of r
        const c = h * h - originalBall.radius * originalBall.radius; // Constant term
        // Use the quadratic formula: r = (-b ± sqrt(b^2 - 4ac)) / 2a
        const discriminant = b * b - 4 * a * c;
        const r1 = (-b + Math.sqrt(discriminant)) / (2 * a);
        const r2 = (-b - Math.sqrt(discriminant)) / (2 * a);
        r = Math.max(r1, r2); // Choose the positive radius

        // Calculate the center point on the plane
        centerPoint = math.add(planePoint, math.multiply(planeNormal, r));
    }

    const offset = math.multiply(planeTangent, r); // Offset in the direction of the tangent

    // Create two smaller balls centered symmetrically about the centerPoint
    const smallerBalls = [
        {
            center: math.add(centerPoint, offset),
            radius: r
        },  // First smaller ball
        {
            center: math.subtract(centerPoint, offset),
            radius: r
        } // Second smaller ball
    ];

    return smallerBalls;
}

function approximateCutBall3(originalBall, planePoint, planeNormal) {
    // Calculate the distance h from the ball center to the plane
    const h = math.dot(math.subtract(originalBall.center, planePoint), planeNormal);

    // Check if the ball intersects the plane
    if (h >= originalBall.radius) {
        return [originalBall]; // Ball is entirely above the plane
    }

    // To find an offset vector parallel to the plane, we can pick one vector in the plane
    let planeTangent1 = math.cross(planeNormal, [1, 0, 0]); // First tangent vector
    let planeTangent2 = math.cross(planeNormal, [0, 1, 0]);
    // Ensure we have a valid tangent, recalculate if the first cross product resulted in a zero vector
    if (math.norm(planeTangent1) === 0) {
        planeTangent1 = math.cross(planeTangent1, [0, 0, 1]);
    }
    if (math.norm(planeTangent2) === 0) {
        planeTangent2 = math.cross(planeTangent2, [0, 0, 1]);
    }

    let r, centerPoint;

    if (h >= math.multiply((2 * math.sqrt(3) + 3) / 3, originalBall.radius)) {
        r = math.multiply((2 * math.sqrt(3) + 3) / 3, originalBall.radius);
        centerPoint = originalBall.center;
    }
    else {
        // Solve for r using the quadratic equation: 4/3 * r^2 + (2R - 2h)r + (h^2 - R^2) = 0
        const a = 4 / 3; // Coefficient of r^2
        const b = 2 * originalBall.radius - 2 * h; // Coefficient of r
        const c = h * h - originalBall.radius * originalBall.radius; // Constant term
        // Use the quadratic formula: r = (-b ± sqrt(b^2 - 4ac)) / 2a
        const discriminant = b * b - 4 * a * c;
        const r1 = (-b + Math.sqrt(discriminant)) / (2 * a);
        const r2 = (-b - Math.sqrt(discriminant)) / (2 * a);
        r = Math.max(r1, r2); // Choose the positive radius

        // Calculate the center point on the plane
        centerPoint = math.add(planePoint, math.multiply(planeNormal, r));
    }

    const offset1 = math.multiply(planeTangent1, 2 / math.sqrt(3) * r); // Offset in the direction of the tangent
    const offset2 = math.add(math.multiply(planeTangent1, -1 / math.sqrt(3) * r), math.multiply(planeTangent2, r));
    const offset3 = math.subtract(math.multiply(planeTangent1, -1 / math.sqrt(3) * r), math.multiply(planeTangent2, r));

    // Create two smaller balls centered symmetrically about the centerPoint
    const smallerBalls = [
        {
            center: math.add(centerPoint, offset1),
            radius: r
        },  // First smaller ball
        {
            center: math.add(centerPoint, offset2),
            radius: r
        }, // Second smaller ball
        {
            center: math.add(centerPoint, offset3),
            radius: r
        } // Second smaller ball
    ];

    return smallerBalls;
}
