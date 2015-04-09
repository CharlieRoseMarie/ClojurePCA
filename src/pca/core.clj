(ns pca.core
  (:gen-class)
  (:use (incanter core stats charts io))
  (:require [clojure.core.reducers :as r]))

(set! *warn-on-reflection* true)

(def THRESH 0.00001)

(defn import-matrix [filename]
  (to-matrix (read-dataset filename :delim \tab)))

(defn normalize [v]
  (let [sum-of (sum-of-squares v)]
    (when-not (zero? sum-of)
      (div v (Math/sqrt sum-of)))))

(defn PCA [m]
  (->> m
       covariance
       decomp-eigenvalue
       sort-components))


(defn sort-components [eigens]
  (let [{:keys [values vectors]} eigens
        coupled-values (map vector values vectors)
        sorted-components (sort-by first > coupled-values)]
    (map second sorted-components)))

(def plastics (import-matrix "resources/plastic_4_10.txt"))
(def p (matrix [[1 2 3 4 5] [2 3 4 5 6] [10 9 8 7 6]]))

(defn calc-component-phase [x u]
  (->> u
       (matrix-projection x)
       normalize
       (matrix-projection (trans x))))

(defn calc-component [m thresh]
  (loop [u (calc-component-phase m (get-nth-col m 0))]
    (if (check-threshold [(get-nth-col m 0) u] thresh)
      u
      (recur (calc-component-phase m u)))))

(defn new-matrix [E t]
  (let [p (normalize (matrix-projection E t))]
    (minus E (mmult t (trans p)))))

(defn fast-PCA [m num-of-components]
  (to-vect
    (rest
      (map first
        (take (inc num-of-components)
               (iterate (fn [[vecs mat]]
                          [(calc-component mat THRESH)
                           (new-matrix mat vecs)])
                        [(get-nth-col m 0) (matrix-mean-center m)]))))))


(calc-component p 0.00001)
(defn check-threshold [values thresh]
  ; Returns false if the difference between the new and old values
  ; is less than the threshold
  (let [[old new] values]
    (> (abs (- (sum-of-squares old)
               (sum-of-squares new)))
       (* (sum-of-squares new) thresh))))

(defn matrix-projection [m v]
  (when-not (zero? (sum-of-squares v))
    (let [projection (mmult (trans m) v)]
      (div projection (sum-of-squares v)))))

(defn get-nth-col [m n]
  (sel m :cols n))

(defn matrix-mean [m]
  (map mean m))

(defn matrix-mean-center [m]
  (let [ma-mean (matrix-mean m)]
    (->> m
       trans
       (map minus ma-mean)
       trans)))
