(ns pca.core
  (:gen-class)
  (:use (incanter core stats charts io))
  (:require [clojure.core.reducers :as r]))

(set! *warn-on-reflection* true)

(def THRESH 0.00001)

(defn import-matrix [filename]
  (to-matrix (read-dataset filename :delim \tab)))

(defn check-threshold [values thresh]
  ; Returns false if the difference between the new and old values
  ; is less than the threshold
  (let [[old new] values]
    (> (abs (- (sum-of-squares old)
               (sum-of-squares new)))
       (* (sum-of-squares new) thresh))))

(defn normalize [v]
  (let [sum-of (sum-of-squares v)]
    (when-not (zero? sum-of)
      (div v (Math/sqrt sum-of)))))

(defn get-nth-col [m n]
  (sel m :cols n))

(defn matrix-projection [m v]
  (when-not (zero? (sum-of-squares v))
    (let [projection (mmult m v)]
      (div projection (sum-of-squares v)))))

(defn calc-scales-and-loadings-pass [m t]
  (let [p (normalize (matrix-projection (trans m) t))]
    {:t (matrix-projection m p) :p p}))

(defn calc-scales-and-loadings [m]
  (let [initial-t (get-nth-col m 0)]
    (loop [sl (calc-scales-and-loadings-pass m initial-t) t initial-t]
      (if (check-threshold [t (:t sl)] THRESH)
        sl
        (recur (calc-scales-and-loadings-pass m {:t sl}) sl)))))

(defn correction-matrix [sl]
  (mmult (:t sl) (trans (:p sl))))

(defn PCA [m n]
  (map first
    (take n
        (iterate (fn [[sl x]]
                   [(calc-scales-and-loadings x) (minus x (correction-matrix sl))])
                   [(calc-scales-and-loadings m) (minus m (correction-matrix (calc-scales-and-loadings m)))]))))

(defn matrix-mean [m]
  (map mean m))

(defn matrix-mean-center [m]
  (let [ma-mean (matrix-mean m)]
    (->> m
       trans
       (map minus ma-mean)
       trans)))
